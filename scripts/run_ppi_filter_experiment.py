"""
scripts/run_ppi_filter_experiment.py — compare RotatE performance under three
levels of BioGRID PPI evidence filtering.

All conditions share the same fixed gene→MP test set (derived from the
unfiltered graph), preventing contamination from independent random splits.

Test-set locking
----------------
1. The canonical G→P test positives are computed once from the *unfiltered*
   graph using the standard 80/10/10 split (seed 42).
2. For each condition, those pairs are pre-assigned to that condition's test
   fold before training, regardless of the condition's entity-to-integer
   mapping (handled by _make_condition_splits in ablation_study.py).
3. All relation types are included in every condition — conditions differ by
   graph content (PPI evidence threshold), not by relation type.

Usage
-----
    uv run python scripts/run_ppi_filter_experiment.py          # full run
    uv run python scripts/run_ppi_filter_experiment.py --fast   # smoke test
    uv run python scripts/run_ppi_filter_experiment.py --out-dir pykeen_out/my_run/
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from src.pykeen.ablation_study import (
    _get_fixed_test_positives,
    _make_condition_splits,
    _run_auprc,
)
from src.pykeen.config import DEFAULT_CONFIG
from src.pykeen.train import run_training
from src.utils.logging_config import setup_logging

logger = logging.getLogger(__name__)

# ── Experiment definitions ────────────────────────────────────────────────────

EXPERIMENTS: dict[str, str] = {
    "ppi_unfiltered": "biocypher_out/pheno_ppi_unfiltered/human/",
    "ppi_filtered_2": "biocypher_out/pheno_ppi_filtered_2/human/",
    "ppi_filtered_3": "biocypher_out/pheno_ppi_filtered_3/human/",
}

# Reference graph for the canonical test split — must be the most complete
# graph so that all G→P pairs are available to the locked test set.
_BASE_DATA_DIR = EXPERIMENTS["ppi_unfiltered"]

_DEFAULT_OUT_DIR = Path("pykeen_out/ppi_filter/")

# Hyperparameters — identical across all conditions so only PPI filtering varies.
# target_steps normalises the gradient budget across conditions regardless of
# graph size (different PPI thresholds produce graphs of different sizes).
_FULL_CONFIG: dict = {
    "target_steps":  50_000,
    "es_patience":   0,     # disabled: full-entity eval every es_frequency epochs
}                           # is prohibitively slow on CPU; target_steps controls budget

# --fast preset: fixed epoch count instead of target_steps so the smoke test
# finishes quickly.  target_steps must be absent — if present, compute_epochs()
# overrides num_epochs and the preset has no effect.
_FAST_CONFIG: dict = {
    "embedding_dim":    128,
    "num_epochs":       100,
    "num_negs_per_pos": 32,
    "batch_size":       2048,
    "es_patience":      0,
}


# ── Main experiment logic ─────────────────────────────────────────────────────

def run_experiment(out_dir: Path, fast: bool = False, model: str = "RotatE") -> pd.DataFrame:
    """Train and evaluate one model per PPI-filter condition.

    Returns a DataFrame with one row per condition.
    Also writes <out_dir>/comparison_summary.csv.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = {**(  _FAST_CONFIG if fast else _FULL_CONFIG), "model": model}

    # ── Fixed test set: derived once from the unfiltered graph ───────────────
    # All conditions are evaluated against exactly the same G→P test positives.
    fixed_positives = _get_fixed_test_positives(_BASE_DATA_DIR)

    rows: list[dict] = []

    for name, data_dir in EXPERIMENTS.items():
        cond_dir = out_dir / name
        logger.info("=" * 60)
        logger.info("Experiment : %s", name)
        logger.info("Data dir   : %s", data_dir)
        logger.info("=" * 60)

        try:
            # Lock canonical test pairs into this condition's test fold.
            # include_relations=None keeps all relation types — conditions differ
            # by graph content (PPI threshold), not by which relations are present.
            train_tf, val_tf, test_tf = _make_condition_splits(
                data_dir=data_dir,
                include_relations=None,
                locked_test_pairs=fixed_positives,
                train_frac=DEFAULT_CONFIG["train_frac"],
                val_frac=DEFAULT_CONFIG["val_frac"],
                seed=DEFAULT_CONFIG["random_seed"],
            )

            summary = run_training(
                data_dir=data_dir,
                out_dir=cond_dir,
                config=cfg,
                train_tf=train_tf,
                val_tf=val_tf,
                test_tf=test_tf,
            )

            auprc_result = _run_auprc(cond_dir, fixed_positives=fixed_positives)

        except Exception:
            logger.exception("Experiment '%s' failed — skipping", name)
            rows.append({"experiment": name, "error": "failed"})
            continue

        m = summary.get("metrics", {})
        rows.append({
            "experiment":    name,
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
    csv_path = out_dir / "comparison_summary.csv"
    df.to_csv(csv_path, index=False)
    logger.info("Comparison summary → %s", csv_path)
    return df


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--out-dir", default=str(_DEFAULT_OUT_DIR),
        help="Root output dir; each condition gets a subdirectory.",
    )
    p.add_argument(
        "--fast", action="store_true",
        help="Smoke-test preset (dim=64, epochs=30). ~10 min per condition on CPU.",
    )
    p.add_argument(
        "--model", default="RotatE",
        choices=["RotatE", "ComplEx", "DistMult"],
        help="KGE model to train (default: RotatE).",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    out_dir = Path(args.out_dir)
    setup_logging(out_dir=out_dir, log_filename="ppi_filter.log")

    df = run_experiment(out_dir=out_dir, fast=args.fast, model=args.model)

    print("\n" + "=" * 70)
    print("PPI FILTER EXPERIMENT SUMMARY")
    print("=" * 70)
    cols = ["experiment", "mean_auprc", "mrr", "hits_at_1", "hits_at_3", "hits_at_10"]
    print(df[[c for c in cols if c in df.columns]].to_string(index=False, float_format="%.4f"))
    print("=" * 70)
