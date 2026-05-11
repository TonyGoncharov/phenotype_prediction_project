"""
run_ablation.py - single entry point for the layer ablation study.

Trains one model per ablation condition, evaluates AUPRC (fixed test set)
and 24-way G→P ranking metrics for each, then writes a comparison table to
<out-dir>/<model>/comparison_summary.csv.

Usage
-----
  # Full run — RotatE, all 5 conditions
  OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 uv run python run_ablation.py

  # Run with ComplEx instead
  uv run python run_ablation.py --model ComplEx --threads 16

  # Smoke test — dim=64, 30 epochs per condition (~10 min/condition on CPU)
  uv run python run_ablation.py --fast --threads 16

  # Single condition
  uv run python run_ablation.py --conditions all pheno_only --threads 16

  # Custom paths
  uv run python run_ablation.py \\
      --data-dir biocypher_out/human/ \\
      --out-dir  pykeen_out/ablation/ \\
      --model RotatE \\
      --threads 16

Thread count
------------
RotatE and ComplEx operate in complex embedding space and are memory-bandwidth
bound on CPU.  The optimal thread count is typically 1–2× physical cores per
socket, not the full logical core count.
Run benchmark_threads.py first to find the sweet spot for your server:

  python src/pykeen/benchmark_threads.py \\
      --data-dir biocypher_out/human/ \\
      --threads 4 8 16 32 64 \\
      --cache-dir .cache/triples \\
      --out benchmark_threads.png
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from src.pykeen.ablation_study import CONDITIONS, run_ablation
from src.utils.logging_config import setup_logging


# ── Pre-flight validation ─────────────────────────────────────────────────────


def check_data_dir(data_dir: Path) -> None:
    """Verify that *data_dir* looks like a valid BioCypher output directory.

    The loader expects *-header.csv + *-part000.csv files produced by
    src/build_graph.py.  If they are absent the ablation will fail deep
    inside training, so we surface the error early with a clear message.
    """
    if not data_dir.exists():
        print(f"ERROR: data directory does not exist: {data_dir.resolve()}")
        print()
        print("Build the human knowledge graph first:")
        print("  uv run python run.py --species human")
        print()
        sys.exit(1)

    headers = list(data_dir.glob("*-header.csv"))
    if not headers:
        print(f"ERROR: no BioCypher CSV files found in {data_dir.resolve()}")
        print()
        print("The directory exists but contains no *-header.csv files.")
        print("Re-run the graph build:")
        print("  uv run python run.py --species human")
        print()
        sys.exit(1)

    print(f"Data directory OK — {len(headers)} relation type(s) found in {data_dir}")
    print()


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
        "--out-dir", default="pykeen_out/ablation/",
        help="Root output dir. Each condition gets its own subdirectory.",
    )
    p.add_argument(
        "--fast", action="store_true",
        help="Smoke-test preset: dim=64, 30 epochs per condition. ~10 min/condition on CPU.",
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
        "--model", default="RotatE",
        help=(
            "PyKEEN model name (e.g. RotatE, ComplEx). "
            "Results go to <out-dir>/<model>/<condition>/. Default: RotatE."
        ),
    )
    p.add_argument(
        "--threads", type=int, default=None,
        help=(
            "Number of CPU threads (sets OMP_NUM_THREADS and related vars). "
            "Must be set before PyTorch is imported, so this flag handles it "
            "for you.  If omitted, PyTorch defaults to os.cpu_count() — "
            "usually suboptimal on large servers.  Run benchmark_threads.py "
            "first to find the sweet spot."
        ),
    )
    return p.parse_args()


# ── Main ──────────────────────────────────────────────────────────────────────


_THREAD_ENV_VARS = [
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "PYTORCH_NUM_THREADS",
]


def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir)
    out_dir  = Path(args.out_dir)

    # Set thread env vars before any PyTorch/OpenMP import happens.
    # torch.set_num_threads() cannot change the OpenMP pool after first import,
    # so the env var is the only reliable way to cap threads.
    if args.threads is not None:
        for var in _THREAD_ENV_VARS:
            os.environ[var] = str(args.threads)
    elif not os.environ.get("OMP_NUM_THREADS"):
        print(
            "WARNING: --threads not set and OMP_NUM_THREADS is unset. "
            f"PyTorch will use all {os.cpu_count()} CPU cores, which is "
            f"usually suboptimal for {args.model} on large servers.\n"
            "Run benchmark_threads.py to find the optimal count, then pass "
            "--threads <N> to this script.\n"
        )

    setup_logging(out_dir=out_dir, log_filename="ablation.log")

    n_threads = args.threads or os.environ.get("OMP_NUM_THREADS") or os.cpu_count()

    print("=" * 70)
    print("LAYER ABLATION STUDY")
    print("=" * 70)
    print(f"Data dir : {data_dir.resolve()}")
    print(f"Out dir  : {out_dir.resolve()}")
    print(f"Model    : {args.model}")
    print(f"Fast mode: {args.fast}")
    print(f"Threads  : {n_threads}")
    print()

    check_data_dir(data_dir)

    selected = (
        {k: CONDITIONS[k] for k in args.conditions}
        if args.conditions else CONDITIONS
    )

    print(f"Conditions to run: {', '.join(selected.keys())}")
    print()

    df = run_ablation(
        data_dir=data_dir,
        out_dir=out_dir,
        conditions=selected,
        fast=args.fast,
        model=args.model,
    )

    # ── Summary table ─────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("ABLATION SUMMARY")
    print("=" * 70)

    primary_cols = [
        "condition", "mean_auprc",
        "gp_mrr", "gp_hits_at_1", "gp_hits_at_3", "gp_hits_at_10",
    ]
    print("\n── Primary metrics (fixed test set / 24-way G→P ranking) ──")
    print(df[[c for c in primary_cols if c in df.columns]].to_string(
        index=False, float_format="%.4f"))

    kge_cols = ["condition", "mrr", "hits_at_1", "hits_at_3", "hits_at_10"]
    print("\n── Standard KGE metrics (entity space differs — not directly comparable) ──")
    print(df[[c for c in kge_cols if c in df.columns]].to_string(
        index=False, float_format="%.4f"))

    print()
    print(f"Full results → {(out_dir / args.model / 'comparison_summary.csv').resolve()}")
    print(f"Log          → {(out_dir / 'ablation.log').resolve()}")
    print("=" * 70)


if __name__ == "__main__":
    main()
