"""
run_ablation.py - single entry point for the layer ablation study.

Trains one RotatE model per ablation condition, evaluates AUPRC against a
fixed test set (no contamination), and writes a comparison table to
<out-dir>/comparison_summary.csv.

The fixed gene→MP test set is computed once from the full graph (seed=42)
inside run_ablation() and locked into every condition's test fold before
training, so AUPRC values are directly comparable across conditions.

Usage
-----
  # Full run — all 5 conditions
  OMP_NUM_THREADS=16 uv run python run_ablation.py

  # Single condition
  uv run python run_ablation.py --conditions all --threads 16

  # Smoke test — dim=64, 30 epochs per condition (~10 min/condition on CPU)
  uv run python run_ablation.py --fast --threads 16

Thread count
------------
RotatE is memory-bandwidth bound on CPU.  The optimal thread count is
typically 1–2× physical cores per socket, not all logical cores.
Run benchmark_threads.py first to find the sweet spot for your server.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from src.pykeen.ablation_study import (
    CONDITIONS,
    run_ablation,
)
from src.utils.logging_config import setup_logging


# ── Pre-flight validation ─────────────────────────────────────────────────────


def check_data_dir(data_dir: Path) -> None:
    if not data_dir.exists():
        print(f"ERROR: data directory does not exist: {data_dir.resolve()}")
        print("\nBuild the human knowledge graph first:")
        print("  uv run python run.py --species human\n")
        sys.exit(1)

    headers = list(data_dir.glob("*-header.csv"))
    if not headers:
        print(f"ERROR: no BioCypher CSV files found in {data_dir.resolve()}")
        print("\nThe directory exists but contains no *-header.csv files.")
        print("Re-run the graph build:")
        print("  uv run python run.py --species human\n")
        sys.exit(1)

    print(f"Data directory OK — {len(headers)} relation type(s) found in {data_dir}\n")


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Layer ablation study for the human gene–phenotype KG.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--data-dir", default="biocypher_out/human/",
        help="BioCypher output directory (produced by run.py --species human).",
    )
    p.add_argument(
        "--out-dir", default=None,
        help="Root output dir. Each condition gets its own subdirectory. "
             "Defaults to pykeen_out/ablation/<model_name_lower>/.",
    )
    p.add_argument(
        "--fast", action="store_true",
        help="Smoke-test preset: dim=64, 30 epochs per condition.",
    )
    p.add_argument(
        "--conditions", nargs="+", default=None,
        choices=list(CONDITIONS.keys()),
        metavar="CONDITION",
        help=(
            "Subset of conditions to run (default: all five). "
            f"Choices: {', '.join(CONDITIONS.keys())}."
        ),
    )
    p.add_argument(
        "--threads", type=int, default=None,
        help=(
            "Number of CPU threads (sets OMP_NUM_THREADS and related vars). "
            "Run benchmark_threads.py first to find the sweet spot."
        ),
    )
    p.add_argument(
        "--model", default="RotatE",
        choices=["RotatE", "ComplEx", "DistMult"],
        help="KGE model to train for every ablation condition.",
    )
    p.add_argument(
        "--target-steps", type=int, default=None,
        help=(
            "Gradient-step budget per condition. Overrides num_epochs so every "
            "condition trains for the same number of parameter updates regardless "
            "of graph size."
        ),
    )
    return p.parse_args()


# ── Thread helpers ────────────────────────────────────────────────────────────

_THREAD_ENV_VARS = [
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "PYTORCH_NUM_THREADS",
]


def _apply_threads(n: int | None) -> None:
    if n is not None:
        for var in _THREAD_ENV_VARS:
            os.environ[var] = str(n)
    elif not os.environ.get("OMP_NUM_THREADS"):
        print(
            f"WARNING: --threads not set and OMP_NUM_THREADS is unset. "
            f"PyTorch will use all {os.cpu_count()} CPU cores, which is "
            "usually suboptimal for KGE models on large servers.\n"
            "Run benchmark_threads.py to find the optimal count.\n"
        )


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> None:
    args     = parse_args()
    data_dir = Path(args.data_dir)
    out_dir  = Path(args.out_dir) if args.out_dir else Path(f"pykeen_out/ablation/{args.model.lower()}")

    _apply_threads(args.threads)
    check_data_dir(data_dir)
    setup_logging(out_dir=out_dir, log_filename="ablation.log")

    n_threads = args.threads or os.environ.get("OMP_NUM_THREADS") or os.cpu_count()

    print("=" * 70)
    print("LAYER ABLATION STUDY")
    print("=" * 70)
    print(f"Data dir  : {data_dir.resolve()}")
    print(f"Out dir   : {out_dir.resolve()}")
    print(f"Fast mode : {args.fast}")
    print(f"Threads   : {n_threads}")
    print()

    selected = (
        {k: CONDITIONS[k] for k in args.conditions}
        if args.conditions else CONDITIONS
    )

    print(f"Conditions to run: {', '.join(selected.keys())}\n")

    df = run_ablation(
        data_dir=data_dir,
        out_dir=out_dir,
        conditions=selected,
        fast=args.fast,
        model=args.model,
        target_steps=args.target_steps,
    )

    # ── Summary table ─────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("ABLATION SUMMARY")
    print("=" * 70)

    cols = ["condition", "mean_auprc", "mrr", "hits_at_1", "hits_at_3", "hits_at_10"]
    print(df[[c for c in cols if c in df.columns]].to_string(index=False, float_format="%.4f"))

    print()
    print(f"Full results → {(out_dir / 'comparison_summary.csv').resolve()}")
    print(f"Log          → {(out_dir / 'ablation.log').resolve()}")
    print("=" * 70)


if __name__ == "__main__":
    main()