from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import torch
from pykeen.pipeline import pipeline
from pykeen.triples import TriplesFactory

from .config import DEFAULT_CONFIG, GENE_PHENOTYPE_RELATION
from .loader import build_triples_factory

logger = logging.getLogger(__name__)


# ── Data splitting ────────────────────────────────────────────────────────────


def split_triples_factory(
    tf: TriplesFactory,
    train_frac: float = DEFAULT_CONFIG["train_frac"],
    val_frac:   float = DEFAULT_CONFIG["val_frac"],
    seed:       int   = DEFAULT_CONFIG["random_seed"],
) -> tuple[TriplesFactory, TriplesFactory, TriplesFactory]:
    """Split *tf* into train / validation / test factories.

    PyKEEN's native split guarantees that every entity and relation that
    appears in the val / test sets also appears in the training set
    (no cold-start entities at evaluation time).

    Returns
    -------
    (train_tf, val_tf, test_tf)
    """
    test_frac = round(1.0 - train_frac - val_frac, 6)
    if test_frac <= 0:
        raise ValueError(
            f"train_frac ({train_frac}) + val_frac ({val_frac}) must be < 1.0"
        )

    train_tf, val_tf, test_tf = tf.split(
        ratios=[train_frac, val_frac, test_frac],
        random_state=seed,
    )

    logger.info(
        "Split — train: %d  |  val: %d  |  test: %d",
        train_tf.num_triples,
        val_tf.num_triples,
        test_tf.num_triples,
    )
    return train_tf, val_tf, test_tf


# ── Training pipeline ─────────────────────────────────────────────────────────


def run_training(
    data_dir:  str | Path,
    out_dir:   str | Path,
    config:    dict | None = None,
    device:    str = "auto",
    cache_dir: str | Path | None = None,
) -> dict:
    """Full RotatE training + evaluation run.

    Parameters
    ----------
    data_dir  : path to biocypher_out/human/ (BioCypher CSV output)
    out_dir   : where to write model, checkpoints, logs, and summary
    config    : dict of hyperparameter overrides (merged with DEFAULT_CONFIG)
    device    : 'cpu', 'cuda', 'mps', or 'auto' (auto-detects GPU)
    cache_dir : directory for caching the TriplesFactory binary.  On the first
                run the factory is built from CSV and saved here; subsequent runs
                load it directly, skipping the CSV parsing step entirely.
                The cache is invalidated automatically when any source CSV is
                newer than the cached file.

    Returns
    -------
    Summary dict with config + evaluation metrics (MRR, Hits@K).
    Also written to <out_dir>/summary.json.
    """
    cfg = {**DEFAULT_CONFIG, **(config or {})}
    include_relations = cfg.pop("include_relations", None)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Load and split data ────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info("Step 1/3 — Loading triples from %s", data_dir)
    logger.info("=" * 60)

    tf = build_triples_factory(
        data_dir,
        include_relations=include_relations,
        cache_dir=cache_dir,
    )
    train_tf, val_tf, test_tf = split_triples_factory(
        tf,
        train_frac=cfg["train_frac"],
        val_frac=cfg["val_frac"],
        seed=cfg["random_seed"],
    )

    # Report how many gene-phenotype triples ended up in each split
    _report_relation_split(tf, train_tf, val_tf, test_tf)

    # ── 2. Resolve device ────────────────────────────────────────────────────
    if device == "auto":
        if torch.cuda.is_available():
            device = "cuda"
        elif torch.backends.mps.is_available():
            device = "mps"
        else:
            device = "cpu"

    # RotatE works in complex space and calls norm() on complex tensors.
    # PyTorch MPS does not yet implement this op → RuntimeError at first batch.
    # Fall back to CPU automatically rather than crashing.
    if device == "mps" and cfg.get("model", "RotatE") == "RotatE":
        logger.warning(
            "MPS (Apple Silicon GPU) does not support complex norm ops "
            "required by RotatE (PyTorch limitation, not a bug in this code). "
            "Falling back to CPU automatically."
        )
        device = "cpu"

    logger.info("Device: %s", device)

    # ── 3. Run PyKEEN pipeline ────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info("Step 2/3 — Training RotatE  (epochs=%d, dim=%d, negs=%d)",
                cfg["num_epochs"], cfg["embedding_dim"], cfg["num_negs_per_pos"])
    logger.info("=" * 60)

    checkpoint_dir = out_dir / "checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)

    result = pipeline(
        # ── Data ─────────────────────────────────────────────────────────────
        training=train_tf,
        validation=val_tf,
        testing=test_tf,

        # ── Model ─────────────────────────────────────────────────────────────
        # RotatE needs even embedding_dim (complex space: dim/2 complex dims).
        model="RotatE",
        model_kwargs=dict(
            embedding_dim=cfg["embedding_dim"],
        ),

        # ── Training loop ─────────────────────────────────────────────────────
        # sLCWA: fast stochastic negative sampling, standard for large KGs.
        training_loop="sLCWA",

        # BasicNegativeSampler: for each positive (h,r,t) draw k random negatives
        # by corrupting either the head or the tail uniformly at random.
        negative_sampler="basic",
        negative_sampler_kwargs=dict(
            num_negs_per_pos=cfg["num_negs_per_pos"],
        ),

        # ── Optimizer ─────────────────────────────────────────────────────────
        optimizer="Adam",
        optimizer_kwargs=dict(lr=cfg["lr"]),

        # ── Epoch / checkpoint config ─────────────────────────────────────────
        training_kwargs=dict(
            num_epochs=cfg["num_epochs"],
            batch_size=cfg["batch_size"],
            checkpoint_name="rotate.pt",
            checkpoint_directory=str(checkpoint_dir),
            checkpoint_frequency=cfg["checkpoint_frequency"],
        ),

        # ── Evaluation ────────────────────────────────────────────────────────
        # Filtered ranking: when ranking a test triple (h,r,t), all other
        # *known* true tails for (h,r,?) are removed from the candidate set.
        # This is the standard protocol for KGE evaluation.
        evaluator="RankBasedEvaluator",
        evaluator_kwargs=dict(filtered=True),
        evaluation_kwargs=dict(batch_size=cfg["eval_batch_size"]),

        # ── Misc ──────────────────────────────────────────────────────────────
        device=device,
        random_seed=cfg["random_seed"],

        # Log per-epoch losses to a CSV file.
        result_tracker="csv",
        result_tracker_kwargs=dict(
            path=str(out_dir / "training_log.csv"),
        ),
    )

    # ── 4. Save artifacts ─────────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info("Step 3/3 — Saving model and results")
    logger.info("=" * 60)

    result.save_to_directory(str(out_dir))
    logger.info("Saved PyKEEN pipeline result to: %s", out_dir)

    # ── 5. Build summary ──────────────────────────────────────────────────────
    summary = {
        "model": cfg["model"],
        "device": device,
        "graph": {
            "num_entities": tf.num_entities,
            "num_relations": tf.num_relations,
            "num_triples_total": tf.num_triples,
            "num_triples_train": train_tf.num_triples,
            "num_triples_val":   val_tf.num_triples,
            "num_triples_test":  test_tf.num_triples,
        },
        "config": cfg,
        "metrics": {
            # 'both' sides, 'realistic' rank type — standard KGE evaluation protocol.
            # Key format: PyKEEN dot-notation (side.rank_type.metric or just metric).
            "hits_at_1":  _safe_metric(result, "hits@1"),
            "hits_at_3":  _safe_metric(result, "hits@3"),
            "hits_at_10": _safe_metric(result, "hits@10"),
            "mrr": _safe_metric(result, "inverse_harmonic_mean_rank"),
            "mr":  _safe_metric(result, "mean_rank"),
        },
    }

    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=str))

    logger.info("─" * 40)
    logger.info("Evaluation results (filtered, realistic):")
    for k, v in summary["metrics"].items():
        logger.info("  %-12s  %s", k, f"{v:.4f}" if isinstance(v, float) else v)
    logger.info("Summary → %s", summary_path)

    return summary


# ── Helpers ───────────────────────────────────────────────────────────────────


def _safe_metric(result, key: str) -> float | None:
    """Extract a scalar metric from a PyKEEN PipelineResult.

    Uses result.get_metric() — the stable API since PyKEEN 1.10.
    Key follows PyKEEN dot-notation, e.g. 'hits@10', 'mean_rank',
    'inverse_harmonic_mean_rank'.  Defaults to both-sides / realistic.
    """
    try:
        return float(result.get_metric(key))
    except Exception:
        return None


def _report_relation_split(
    full_tf:    TriplesFactory,
    train_tf:   TriplesFactory,
    val_tf:     TriplesFactory,
    test_tf:    TriplesFactory,
) -> None:
    """Log how many gene→phenotype triples appear in each split."""
    if GENE_PHENOTYPE_RELATION not in full_tf.relation_to_id:
        return
    rel_id = full_tf.relation_to_id[GENE_PHENOTYPE_RELATION]

    def _count(tf: TriplesFactory) -> int:
        # mapped_triples is a (N,3) int tensor: [head_id, relation_id, tail_id]
        return int((tf.mapped_triples[:, 1] == rel_id).sum().item())

    logger.info(
        "Gene→Phenotype triples — train: %d  val: %d  test: %d",
        _count(train_tf), _count(val_tf), _count(test_tf),
    )


# ── CLI entry point ───────────────────────────────────────────────────────────


_FAST_CONFIG = {
    "embedding_dim":    64,
    "num_epochs":       30,
    "num_negs_per_pos": 32,
    "batch_size":       1024,
}


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Train RotatE on the Human Knowledge Graph (gene→phenotype focus).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--data-dir", required=True,
        help="BioCypher output dir that contains *-header.csv and *-part000.csv files "
             "(e.g. biocypher_out/human/).",
    )
    p.add_argument(
        "--out-dir", default="pykeen_out/rotate",
        help="Directory to save the trained model, logs, and summary.",
    )
    p.add_argument("--epochs",  type=int,   default=DEFAULT_CONFIG["num_epochs"])
    p.add_argument("--dim",     type=int,   default=DEFAULT_CONFIG["embedding_dim"],
                   help="Embedding dimension (must be even for RotatE).")
    p.add_argument("--lr",      type=float, default=DEFAULT_CONFIG["lr"])
    p.add_argument("--batch",   type=int,   default=DEFAULT_CONFIG["batch_size"])
    p.add_argument("--negs",    type=int,   default=DEFAULT_CONFIG["num_negs_per_pos"],
                   help="Negative samples per positive triple.")
    p.add_argument("--device",  default="auto",
                   choices=["auto", "cpu", "cuda", "mps"])
    p.add_argument("--seed",    type=int,   default=DEFAULT_CONFIG["random_seed"])
    p.add_argument(
        "--fast", action="store_true",
        help="Smoke-test preset: dim=64, epochs=20, negs=32, batch=1024. "
             "Completes in ~10 min on CPU. Verify the pipeline before a full run.",
    )
    p.add_argument(
        "--cache-dir", default=None,
        help="Directory for caching the TriplesFactory binary. "
             "Skips CSV parsing on subsequent runs with the same data. "
             "Cache is invalidated automatically when source CSVs change.",
    )
    return p.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )
    args = _parse_args()

    if args.fast:
        override = _FAST_CONFIG.copy()
        logger.info("--fast preset: %s", override)
    else:
        override = {
            "num_epochs":       args.epochs,
            "embedding_dim":    args.dim,
            "lr":               args.lr,
            "batch_size":       args.batch,
            "num_negs_per_pos": args.negs,
            "random_seed":      args.seed,
        }

    dim = override.get("embedding_dim", args.dim)
    if dim % 2 != 0:
        raise ValueError(f"embedding_dim must be even for RotatE (got {dim})")

    run_training(
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        config=override,
        device=args.device,
        cache_dir=args.cache_dir,
    )