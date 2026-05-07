"""benchmark_threads.py — find the optimal OMP_NUM_THREADS for RotatE on this server.

Each thread count is tested in a **separate subprocess** so that OMP_NUM_THREADS
is set before any PyTorch / OpenMP import.  This is the only correct way:
torch.set_num_threads() cannot change the OpenMP thread pool after first import.

Usage
-----
    python benchmark_threads.py \\
        --data-dir biocypher_out/human/ \\
        --out benchmark_threads.png

    # Test only specific counts (faster):
    python benchmark_threads.py \\
        --data-dir biocypher_out/human/ \\
        --threads 8 16 32 64 128 \\
        --out benchmark_threads.png
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

logger = logging.getLogger(__name__)

# ── Benchmark training config ─────────────────────────────────────────────────
# 5 epochs is enough to measure throughput; small enough to finish in minutes.
_BENCH_EPOCHS   = 5
_BENCH_DIM      = 64
_BENCH_BATCH    = 1024
_BENCH_NEGS     = 32


# ── Thread count helpers ──────────────────────────────────────────────────────

def _default_thread_counts() -> list[int]:
    """Powers of 2 from 1 up to os.cpu_count(), always including cpu_count."""
    max_n = os.cpu_count() or 1
    counts: list[int] = []
    n = 1
    while n <= max_n:
        counts.append(n)
        n *= 2
    if counts[-1] != max_n:
        counts.append(max_n)
    return counts


# ── Single subprocess run ─────────────────────────────────────────────────────

def _run_one(
    n_threads: int,
    data_dir: Path,
    out_dir: Path,
    cache_dir: Path | None,
) -> float | None:
    """Launch one training run in a subprocess with OMP_NUM_THREADS=n_threads.

    Returns wall-clock seconds, or None on failure.
    """
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(n_threads)
    env["MKL_NUM_THREADS"] = str(n_threads)

    cmd = [
        sys.executable, "-m", "src.pykeen.train",
        "--data-dir", str(data_dir),
        "--out-dir",  str(out_dir),
        "--epochs",   str(_BENCH_EPOCHS),
        "--dim",      str(_BENCH_DIM),
        "--batch",    str(_BENCH_BATCH),
        "--negs",     str(_BENCH_NEGS),
        "--device",   "cpu",
    ]
    if cache_dir is not None:
        cmd += ["--cache-dir", str(cache_dir)]

    t0 = time.perf_counter()
    try:
        result = subprocess.run(
            cmd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=3600,
        )
    except subprocess.TimeoutExpired:
        logger.error("Run with threads=%d timed out after 3600s", n_threads)
        return None
    elapsed = time.perf_counter() - t0

    if result.returncode != 0:
        logger.error("Run with threads=%d failed:\n%s", n_threads, result.stdout[-2000:])
        return None

    return elapsed


# ── Benchmark loop ────────────────────────────────────────────────────────────

def run_benchmark(
    data_dir: Path,
    thread_counts: list[int],
    tmp_root: Path,
    cache_dir: Path | None = None,
) -> dict[int, float]:
    """Run one subprocess per thread count; return {n_threads: wall_seconds}."""
    results: dict[int, float] = {}

    print(f"CPU count : {os.cpu_count()}")
    print(f"Threads   : {thread_counts}")
    print(f"Epochs    : {_BENCH_EPOCHS}  dim={_BENCH_DIM}  batch={_BENCH_BATCH}")
    print(f"Cache dir : {cache_dir or 'disabled (slow — each run re-parses CSV)'}")
    print()

    for n in thread_counts:
        out_dir = tmp_root / f"bench_{n}"
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"  threads={n:4d} ...", end="", flush=True)
        elapsed = _run_one(n, data_dir, out_dir, cache_dir)

        if elapsed is None:
            print("  FAILED")
        else:
            results[n] = elapsed
            print(f"  {elapsed:6.1f}s")

    return results


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_results(results: dict[int, float], out_path: Path) -> None:
    if not results:
        print("No results to plot.")
        return

    threads  = sorted(results)
    times    = [results[n] for n in threads]
    best_n   = min(results, key=results.__getitem__)
    baseline = results[threads[0]]
    speedups = [baseline / results[n] for n in threads]
    ideal    = [min(n / threads[0], max(speedups)) for n in threads]

    fig, (ax_t, ax_s) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"OMP_NUM_THREADS benchmark — RotatE  "
        f"({_BENCH_EPOCHS} epochs, dim={_BENCH_DIM}, batch={_BENCH_BATCH})",
        fontsize=13, fontweight="bold",
    )

    # ── Left: wall time ───────────────────────────────────────────────────── #
    colors = ["#2ecc71" if n == best_n else "#4a90d9" for n in threads]
    bars = ax_t.bar([str(n) for n in threads], times, color=colors, width=0.6)
    for bar, t in zip(bars, times):
        ax_t.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(times) * 0.01,
            f"{t:.1f}s", ha="center", va="bottom", fontsize=9,
        )
    ax_t.set_xlabel("OMP_NUM_THREADS", fontsize=11)
    ax_t.set_ylabel("Wall time (seconds)", fontsize=11)
    ax_t.set_title("Training time (lower = faster)")
    ax_t.set_ylim(0, max(times) * 1.18)
    ax_t.legend(handles=[
        Patch(facecolor="#2ecc71", label=f"Best: {best_n} threads"),
        Patch(facecolor="#4a90d9", label="Other"),
    ], fontsize=9)
    ax_t.grid(axis="y", alpha=0.25)

    # ── Right: speedup ────────────────────────────────────────────────────── #
    xs = [str(n) for n in threads]
    ax_s.plot(xs, speedups, "o-", color="#4a90d9", lw=2, ms=7, label="Actual")
    ax_s.plot(xs, ideal,    "--", color="#bbb",    lw=1.5,      label="Ideal (linear)")
    ax_s.axvline(
        x=str(best_n), color="#2ecc71",
        linestyle=":", lw=2, alpha=0.9, label=f"Best: {best_n} threads",
    )
    for i, (n, s) in enumerate(zip(threads, speedups)):
        ax_s.annotate(
            f"{s:.2f}×", xy=(i, s),
            xytext=(0, 8), textcoords="offset points",
            ha="center", fontsize=8,
        )
    ax_s.set_xlabel("OMP_NUM_THREADS", fontsize=11)
    ax_s.set_ylabel(f"Speedup vs {threads[0]} thread(s)", fontsize=11)
    ax_s.set_title("Speedup")
    ax_s.set_ylim(bottom=0)
    ax_s.legend(fontsize=9)
    ax_s.grid(axis="y", alpha=0.25)

    plt.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")

    # ── Console summary ───────────────────────────────────────────────────── #
    print("\nResults:")
    for n in threads:
        marker = "  ← best" if n == best_n else ""
        print(f"  {n:4d} threads:  {results[n]:.1f}s{marker}")

    print(f"\nPlot saved: {out_path}")
    print(f"\nRecommendation:")
    print(f"  OMP_NUM_THREADS={best_n} MKL_NUM_THREADS={best_n} \\")
    print(f"  uv run python -m src.pykeen.train --data-dir biocypher_out/human/ ...")

    # Save raw numbers alongside the chart
    json_path = out_path.with_suffix(".json")
    json_path.write_text(json.dumps(
        {"thread_counts": threads, "wall_seconds": times, "best_n_threads": best_n},
        indent=2,
    ))
    print(f"Raw data  : {json_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Benchmark OMP_NUM_THREADS for RotatE (subprocess-based).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--data-dir", required=True,
                   help="BioCypher output dir, e.g. biocypher_out/human/")
    p.add_argument("--out", default="benchmark_threads.png",
                   help="Output PNG path.")
    p.add_argument("--threads", type=int, nargs="+", default=None,
                   help="Explicit thread counts to test. "
                        "Default: powers of 2 up to os.cpu_count().")
    p.add_argument("--cache-dir", default=".cache/triples",
                   help="TriplesFactory cache dir. Strongly recommended: without it "
                        "each subprocess re-parses all CSV files (~1-2 min each).")
    p.add_argument("--tmp-dir", default=None,
                   help="Temp dir for per-run outputs (auto-created if omitted).")
    return p.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )

    args  = _parse_args()
    counts = args.threads or _default_thread_counts()

    own_tmp = args.tmp_dir is None
    tmp_root = Path(args.tmp_dir) if args.tmp_dir else Path(tempfile.mkdtemp(prefix="bench_threads_"))

    try:
        results = run_benchmark(
            data_dir=Path(args.data_dir),
            thread_counts=counts,
            tmp_root=tmp_root,
            cache_dir=Path(args.cache_dir) if args.cache_dir else None,
        )
        plot_results(results, out_path=Path(args.out))
    finally:
        if own_tmp and tmp_root.exists():
            shutil.rmtree(tmp_root, ignore_errors=True)