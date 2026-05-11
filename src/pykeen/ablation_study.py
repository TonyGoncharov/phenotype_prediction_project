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
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

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

_PHENO = {"HumanGeneHasMpTopTerm"}
_GO    = {"HumanGeneHasGoTerm"}
# TissueMappedToUberon ties tissue nodes to the Uberon hierarchy;
# grouped with expression since both are absent when GTEx data is excluded.
_EXPR  = {"HumanGeneExpressedInTissue", "TissueMappedToUberon"}
_PPI   = {"HumanGeneEncodesProtein", "HumanProteinInteractsWith"}

# ── Ablation conditions ───────────────────────────────────────────────────────
# include_relations=None  →  all relations present in the graph.
# An explicit set restricts training to those relation types only.

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
    data_dir:   str | Path = "biocypher_out/human/",
    out_dir:    str | Path = "pykeen_out/ablation/",
    conditions: dict[str, set[str] | None] | None = None,
    fast:       bool = False,
    model:      str  = "RotatE",
) -> pd.DataFrame:
    """Train and evaluate one model per ablation condition.

    Returns a DataFrame with one row per condition (suitable for plotting).
    Also writes <out_dir>/<model>/comparison_summary.csv.

    All conditions share the same fixed gene→phenotype test set (derived from
    the full-graph split) so that AUPRC values are directly comparable.

    Note on MRR/Hits@K: these metrics rank against *all entities in each
    condition's graph*, so the candidate pool differs across conditions
    (e.g. removing PPI shrinks entities from 60k → 36k, artificially boosting
    ranking metrics).  MRR/Hits@K are reported for completeness but should not
    be used as the primary cross-condition comparison metric.
    """
    conditions = conditions or CONDITIONS
    out_path   = Path(out_dir)
    model_path = out_path / model
    model_path.mkdir(parents=True, exist_ok=True)

    # ── Compute fixed test positives ONCE from the full graph ─────────────────
    fixed_positives = _get_fixed_test_positives(data_dir)

    rows: list[dict] = []

    for name, rels in conditions.items():
        cond_dir  = model_path / name
        rels_desc = "all" if rels is None else "|".join(sorted(rels))

        logger.info("=" * 60)
        logger.info("Condition : %s", name)
        logger.info("Relations : %s", rels_desc)
        logger.info("=" * 60)

        try:
            cfg = {**(fast and _FAST_OVERRIDE or {}), "include_relations": rels, "model": model}
            summary = run_training(data_dir=str(data_dir), out_dir=cond_dir, config=cfg)
            _patch_summary(cond_dir, summary)

            auprc_result = _run_auprc(cond_dir, fixed_positives=fixed_positives)

        except Exception:
            logger.exception("Condition '%s' failed — skipping", name)
            rows.append({"condition": name, "error": "failed"})
            continue

        m = summary.get("metrics", {})
        gp = summary.get("gp_metrics", {})
        rows.append({
            "condition":     name,
            "relations":     rels_desc,
            "mean_auprc":    auprc_result.get("mean_auprc"),
            "auprc_classes": auprc_result.get("num_classes"),
            "mrr":           m.get("mrr"),
            "hits_at_1":     m.get("hits_at_1"),
            "hits_at_3":     m.get("hits_at_3"),
            "hits_at_10":    m.get("hits_at_10"),
            # 24-way tail-only metrics — directly comparable across conditions
            "gp_mrr":        gp.get("gp_mrr"),
            "gp_hits_at_1":  gp.get("gp_hits_at_1"),
            "gp_hits_at_3":  gp.get("gp_hits_at_3"),
            "gp_hits_at_10": gp.get("gp_hits_at_10"),
            "gp_n_candidates": gp.get("gp_n_candidates"),
        })
        logger.info(
            "Done — AUPRC=%.4f  MRR=%.4f  Hits@10=%.4f",
            auprc_result.get("mean_auprc") or 0,
            m.get("mrr") or 0,
            m.get("hits_at_10") or 0,
        )

    df = pd.DataFrame(rows)
    csv_path = model_path / "comparison_summary.csv"
    df.to_csv(csv_path, index=False)
    logger.info("Comparison summary → %s", csv_path)
    return df


# ── Helpers ───────────────────────────────────────────────────────────────────

def _run_auprc(
    cond_dir: Path,
    fixed_positives: "dict[str, set[str]] | None" = None,
) -> dict:
    """Load test_triples, compute AUPRC, write per-class CSV, patch summary.json.

    If *fixed_positives* is provided it is forwarded to evaluate_auprc so that
    all conditions are scored against the same canonical test set.
    """
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
    """Write summary dict to <cond_dir>/summary.json (run_training already does this,
    but we overwrite to ensure include_relations is recorded)."""
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
                   help="Root output dir. Results go to <out-dir>/<model>/<condition>/.")
    p.add_argument("--model", default="RotatE",
                   help="PyKEEN model name (e.g. RotatE, ComplEx). Default: RotatE.")
    p.add_argument("--fast", action="store_true",
                   help="Smoke-test preset (dim=64, epochs=30). ~10 min per condition on CPU.")
    p.add_argument("--conditions", nargs="+", default=None,
                   choices=list(CONDITIONS.keys()),
                   help="Subset of conditions to run (default: all).")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()

    # Console + file logging; log file goes next to comparison_summary.csv.
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
        model=args.model,
    )

    print("\n" + "=" * 70)
    print("ABLATION SUMMARY")
    print("=" * 70)

    # Primary metrics (AUPRC fixed test set; gp_* 24-way comparable ranking)
    primary_cols = [
        "condition", "mean_auprc",
        "gp_mrr", "gp_hits_at_1", "gp_hits_at_3", "gp_hits_at_10",
    ]
    print("\n── Primary metrics (fixed test set / 24-way) ──")
    print(df[[c for c in primary_cols if c in df.columns]].to_string(
        index=False, float_format="%.4f"))

    # Standard KGE metrics (informational; NOT comparable across conditions)
    kge_cols = ["condition", "mrr", "hits_at_1", "hits_at_3", "hits_at_10"]
    print("\n── Standard KGE metrics (entity space differs — not directly comparable) ──")
    print(df[[c for c in kge_cols if c in df.columns]].to_string(
        index=False, float_format="%.4f"))

    print("=" * 70)