"""src/pykeen/ablation_study.py — layer ablation study for the human KG.

Trains one RotatE model per condition, evaluates AUPRC + MRR/Hits@K for each,
and writes a comparison table to <out-dir>/comparison_summary.csv.

Usage
-----
    # Full training (~hours per condition on CPU)
    uv run python -m src.pykeen.ablation_study

    # Smoke test (30 epochs, dim=64)
    uv run python -m src.pykeen.ablation_study --fast

    # Single condition
    uv run python -m src.pykeen.ablation_study --conditions all pheno_only

Test-set integrity
------------------
A fixed gene→MP test set is computed once from the full graph (seed=42).
Each condition's TriplesFactory is built by:
  1. Loading the condition's triples.
  2. Pre-assigning the fixed test pairs to the test fold.
  3. Splitting the remainder into train / val at the original 80:10 ratio.

This prevents the contamination that arises when every condition does an
independent random split: because entity→integer mappings differ across
conditions, the same seed produces different folds, and some test pairs
leak into a condition's training set (worst case: pheno_only).
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import torch
import pandas as pd
from pykeen.triples import TriplesFactory

from ..utils.logging_config import setup_logging
from .config import DEFAULT_CONFIG, GENE_PHENOTYPE_RELATION
from .evaluate import evaluate_auprc
from .loader import build_triples_factory
from .predict import GenePhenotypePredictor
from .train import run_training, split_triples_factory

logger = logging.getLogger(__name__)

# ── Relation groups ───────────────────────────────────────────────────────────

_PHENO = {GENE_PHENOTYPE_RELATION}           # was: {"HumanGeneHasMpTopTerm"}
_GO    = {"HumanGeneHasGoTerm"}
_EXPR  = {"HumanGeneExpressedInTissue", "TissueMappedToUberon"}
_PPI   = {"HumanGeneEncodesProtein", "HumanProteinInteractsWith"}

# ── Ablation conditions ───────────────────────────────────────────────────────

CONDITIONS: dict[str, set[str] | None] = {
    "all":        None,
    "no_go":      _PHENO | _EXPR | _PPI,
    "no_ppi":     _PHENO | _GO   | _EXPR,
    "no_expr":    _PHENO | _GO   | _PPI,
    "pheno_only": _PHENO,
}

_FAST_OVERRIDE = {
    "embedding_dim":    64,
    "num_epochs":       30,
    "num_negs_per_pos": 32,
    "batch_size":       1024,
}


# ── Contamination-safe split ──────────────────────────────────────────────────

def _make_condition_splits(
    data_dir:          str | Path,
    include_relations: set[str] | None,
    locked_test_pairs: dict[str, set[str]],
    train_frac:        float = DEFAULT_CONFIG["train_frac"],
    val_frac:          float = DEFAULT_CONFIG["val_frac"],
    seed:              int   = DEFAULT_CONFIG["random_seed"],
) -> tuple[TriplesFactory, TriplesFactory, TriplesFactory]:
    """Return (train_tf, val_tf, test_tf) with zero contamination.

    Gene→MP pairs that appear in *locked_test_pairs* are pre-assigned to the
    test fold and never seen during training, regardless of the condition's
    entity-to-integer mapping.

    Parameters
    ----------
    locked_test_pairs:
        mp_term_id → {gene_id, …}  — the canonical test positives from the
        full-graph split (output of ``_get_fixed_test_positives``).
    """
    tf = build_triples_factory(data_dir, include_relations=include_relations)
    rel_id = tf.relation_to_id.get(GENE_PHENOTYPE_RELATION)

    if rel_id is None:
        # This condition has no G→P relation at all (shouldn't happen in
        # practice, but be safe).  Fall through to a standard split so
        # training can still run.
        logger.warning(
            "'%s' absent from condition — falling back to standard split",
            GENE_PHENOTYPE_RELATION,
        )
        return split_triples_factory(tf, train_frac, val_frac, seed)

    id_to_entity: dict[int, str] = {v: k for k, v in tf.entity_to_id.items()}

    # ── Find which triples in this TriplesFactory are locked test pairs ───────
    locked_mask = torch.zeros(tf.num_triples, dtype=torch.bool)
    for i, (h, r, t) in enumerate(tf.mapped_triples.tolist()):
        if r == rel_id:
            gene   = id_to_entity.get(h, "")
            mp_trm = id_to_entity.get(t, "")
            if mp_trm in locked_test_pairs and gene in locked_test_pairs[mp_trm]:
                locked_mask[i] = True

    n_locked = int(locked_mask.sum().item())
    logger.info(
        "Condition split: %d / %d G→P triples locked into test fold",
        n_locked, int((tf.mapped_triples[:, 1] == rel_id).sum().item()),
    )

    if n_locked == 0:
        # None of the fixed test pairs are present in this condition (e.g.
        # all G→P edges were filtered).  Standard split is safe.
        logger.warning(
            "No locked test pairs found in this condition — standard split used."
        )
        return split_triples_factory(tf, train_frac, val_frac, seed)

    # ── Partition: locked → test, remainder → train + val ────────────────────
    locked_triples    = tf.mapped_triples[locked_mask]
    remainder_triples = tf.mapped_triples[~locked_mask]

    # Reconstruct TriplesFactory for the remainder, preserving the shared
    # entity/relation vocabulary so embeddings are compatible across conditions.
    remainder_tf = TriplesFactory(
        mapped_triples=remainder_triples,
        entity_to_id=tf.entity_to_id,
        relation_to_id=tf.relation_to_id,
    )

    # Split remainder into train / val at the original ratio.
    # Pass a single ratio < 1.0 so PyKEEN returns exactly two splits.
    adj_train_ratio = train_frac / (train_frac + val_frac)   # e.g. 0.8/0.9 ≈ 0.889
    train_tf, val_tf = remainder_tf.split(
        ratios=[adj_train_ratio],
        random_state=seed,
    )

    # Test factory from locked triples — same shared vocabulary.
    test_tf = TriplesFactory(
        mapped_triples=locked_triples,
        entity_to_id=tf.entity_to_id,
        relation_to_id=tf.relation_to_id,
    )

    logger.info(
        "Contamination-safe split — train: %d  val: %d  test: %d",
        train_tf.num_triples, val_tf.num_triples, test_tf.num_triples,
    )
    return train_tf, val_tf, test_tf


# ── Core logic ────────────────────────────────────────────────────────────────

def _get_fixed_test_positives(
    data_dir:   str | Path,
    train_frac: float = DEFAULT_CONFIG["train_frac"],
    val_frac:   float = DEFAULT_CONFIG["val_frac"],
    seed:       int   = DEFAULT_CONFIG["random_seed"],
) -> dict[str, set[str]]:
    """Build the canonical gene→phenotype test set from the *full* graph.

    Loads all relation types, splits with the given seed, and returns the
    gene→phenotype pairs that landed in the test fold as a dict:
        { mp_term_id (str) → set of gene_id strings }

    This is computed once before the ablation loop so that every condition is
    evaluated against exactly the same positive pairs, making AUPRC scores
    directly comparable regardless of which auxiliary edge types are present.
    """
    logger.info("Computing fixed G→P test positives from full graph …")
    tf = build_triples_factory(data_dir, include_relations=None)
    _, _, test_tf = split_triples_factory(
        tf, train_frac=train_frac, val_frac=val_frac, seed=seed
    )

    rel_id = test_tf.relation_to_id.get(GENE_PHENOTYPE_RELATION)
    if rel_id is None:
        raise ValueError(
            f"'{GENE_PHENOTYPE_RELATION}' not found in the full graph. "
            "Cannot build fixed test positives."
        )

    mask = test_tf.mapped_triples[:, 1] == rel_id
    id_to_entity = {v: k for k, v in test_tf.entity_to_id.items()}

    positives: dict[str, set[str]] = {}
    for head_id, _, tail_id in test_tf.mapped_triples[mask].tolist():
        gene    = id_to_entity[head_id]
        mp_term = id_to_entity[tail_id]
        positives.setdefault(mp_term, set()).add(gene)

    total = sum(len(v) for v in positives.values())
    logger.info(
        "Fixed test set: %d gene→MP pairs across %d MP terms",
        total, len(positives),
    )
    return positives


def run_ablation(
    data_dir: str | Path = "biocypher_out/human/",
    out_dir:  str | Path = "pykeen_out/ablation/",
    conditions: dict[str, set[str] | None] | None = None,
    fast: bool = False,
) -> pd.DataFrame:
    """Train and evaluate one model per ablation condition.

    Returns a DataFrame with one row per condition (suitable for plotting).
    Also writes <out_dir>/comparison_summary.csv.

    All conditions share the same fixed gene→phenotype test set (derived from
    the full-graph split) so that AUPRC values are directly comparable.
    Contamination is prevented by pre-assigning the locked test pairs to each
    condition's test fold before training (see ``_make_condition_splits``).

    Note on MRR/Hits@K: these metrics rank against *all entities in each
    condition's graph*, so the candidate pool differs across conditions
    (e.g. removing PPI shrinks entities from 60k → 36k, artificially boosting
    ranking metrics).  MRR/Hits@K are reported for completeness but should not
    be used as the primary cross-condition comparison metric.
    """
    conditions = conditions or CONDITIONS
    out_path   = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # ── Compute fixed test positives ONCE from the full graph ─────────────────
    fixed_positives = _get_fixed_test_positives(data_dir)

    rows: list[dict] = []

    for name, rels in conditions.items():
        cond_dir  = out_path / name
        rels_desc = "all" if rels is None else "|".join(sorted(rels))

        logger.info("=" * 60)
        logger.info("Condition : %s", name)
        logger.info("Relations : %s", rels_desc)
        logger.info("=" * 60)

        try:
            cfg = {**(fast and _FAST_OVERRIDE or {}), "include_relations": rels}

            # ── Contamination-safe splits ─────────────────────────────────────
            # Each condition's test fold is locked to exactly the full-graph
            # test pairs; the remainder is split train/val at the adjusted ratio.
            train_tf, val_tf, test_tf = _make_condition_splits(
                data_dir=data_dir,
                include_relations=rels,
                locked_test_pairs=fixed_positives,
                train_frac=DEFAULT_CONFIG["train_frac"],
                val_frac=DEFAULT_CONFIG["val_frac"],
                seed=DEFAULT_CONFIG["random_seed"],
            )

            summary = run_training(
                data_dir=str(data_dir),
                out_dir=cond_dir,
                config=cfg,
                train_tf=train_tf,
                val_tf=val_tf,
                test_tf=test_tf,
            )
            _patch_summary(cond_dir, summary)

            auprc_result = _run_auprc(cond_dir, fixed_positives=fixed_positives)

        except Exception:
            logger.exception("Condition '%s' failed — skipping", name)
            rows.append({"condition": name, "error": "failed"})
            continue

        m = summary.get("metrics", {})
        rows.append({
            "condition":     name,
            "relations":     rels_desc,
            "mean_auprc":    auprc_result.get("mean_auprc"),
            "auprc_classes": auprc_result.get("num_classes"),
            "mrr":           m.get("mrr"),
            "hits_at_1":     m.get("hits_at_1"),
            "hits_at_3":     m.get("hits_at_3"),
            "hits_at_10":    m.get("hits_at_10"),
        })
        logger.info(
            "Done — AUPRC=%.4f  MRR=%.4f  Hits@10=%.4f",
            auprc_result.get("mean_auprc") or 0,
            m.get("mrr") or 0,
            m.get("hits_at_10") or 0,
        )

    df = pd.DataFrame(rows)
    csv_path = out_path / "comparison_summary.csv"
    df.to_csv(csv_path, index=False)
    logger.info("Comparison summary → %s", csv_path)
    return df


# ── Helpers ───────────────────────────────────────────────────────────────────

def _run_auprc(
    cond_dir: Path,
    fixed_positives: "dict[str, set[str]] | None" = None,
) -> dict:
    """Load test_triples, compute AUPRC, write per-class CSV, patch summary.json."""
    test_tf_path = cond_dir / "test_triples"
    if not test_tf_path.exists():
        logger.warning("test_triples/ not found in %s — AUPRC skipped", cond_dir)
        return {"mean_auprc": None, "num_classes": 0}

    predictor = GenePhenotypePredictor.from_directory(cond_dir)
    test_tf   = TriplesFactory.from_path_binary(test_tf_path)
    result    = evaluate_auprc(predictor, test_tf, override_positives=fixed_positives)

    result["per_class"].to_csv(cond_dir / "auprc_per_phenotype.csv", index=False)

    sp = cond_dir / "summary.json"
    if sp.exists():
        s = json.loads(sp.read_text())
        s["metrics"]["mean_auprc"]    = result["mean_auprc"]
        s["metrics"]["auprc_classes"] = result["num_classes"]
        sp.write_text(json.dumps(s, indent=2, default=str))

    return result


def _patch_summary(cond_dir: Path, summary: dict) -> None:
    sp = cond_dir / "summary.json"
    if sp.exists():
        sp.write_text(json.dumps(summary, indent=2, default=str))


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Layer ablation study for the human gene–phenotype KG.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--data-dir", default="biocypher_out/human/",
                   help="BioCypher output directory.")
    p.add_argument("--out-dir",  default="pykeen_out/ablation/",
                   help="Root output dir; each condition gets a subdirectory.")
    p.add_argument("--fast", action="store_true",
                   help="Smoke-test preset (dim=64, epochs=30). ~10 min per condition on CPU.")
    p.add_argument("--conditions", nargs="+", default=None,
                   choices=list(CONDITIONS.keys()),
                   help="Subset of conditions to run (default: all).")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()

    setup_logging(
        out_dir=Path(args.out_dir),
        log_filename="ablation.log",
    )
    selected = (
        {k: CONDITIONS[k] for k in args.conditions}
        if args.conditions else CONDITIONS
    )

    df = run_ablation(
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        conditions=selected,
        fast=args.fast,
    )

    print("\n" + "=" * 70)
    print("ABLATION SUMMARY")
    print("=" * 70)
    cols = ["condition", "mean_auprc", "mrr", "hits_at_1", "hits_at_3", "hits_at_10"]
    print(df[[c for c in cols if c in df.columns]].to_string(index=False, float_format="%.4f"))
    print("=" * 70)