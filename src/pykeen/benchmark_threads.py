"""benchmark_threads.py — find the optimal OMP_NUM_THREADS for RotatE on this server.

Runs a short fixed-epoch training loop for each thread count, measures wall-clock
time, and saves a two-panel PNG chart.

Usage
-----
    uv run python -m src.pykeen.benchmark_threads \\
        --data-dir biocypher_out/ \\
        --cache-dir .cache/triples \\
        --out benchmark_threads.png

Note: torch.set_num_threads() controls PyTorch's intraop thread pool at runtime.
For OpenMP-backed ops, OMP_NUM_THREADS must be set before the first torch import.
To benchmark that layer too, wrap this script:
    for T in 1 4 8 16 32 64; do
        OMP_NUM_THREADS=$T uv run python -m src.pykeen.benchmark_threads ...
    done
"""

from __future__ import annotations

import argparse
import logging
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no display needed on a server
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import torch

from src.pykeen.train import run_training

# ── Benchmark config ──────────────────────────────────────────────────────────
# Enough epochs to get a stable time signal; small enough to finish quickly.
_BENCH_CONFIG = {
    "num_epochs":        5,
    "embedding_dim":     64,
    "batch_size":        1024,
    "num_negs_per_pos":  32,
    "checkpoint_frequency": 0,   # no checkpoints — pure training time
    "eval_batch_size":   512,
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def _thread_counts(max_n: int) -> list[int]:
    """Powers of 2 from 1 up to max_n, always including max_n."""
    counts: list[int] = []
    n = 1
    while n <= max_n:
        counts.append(n)
        n *= 2
    if counts[-1] != max_n:
        counts.append(max_n)
    return counts


def _silence_training_logs() -> None:
    """Reduce log noise during benchmark runs — keep only WARNING+."""
    for name in ("src.pykeen", "pykeen", "kg_pipeline"):
        logging.getLogger(name).setLevel(logging.WARNING)


# ── Benchmark runner ──────────────────────────────────────────────────────────

def run_benchmark(
    data_dir: str | Path,
    cache_dir: str | Path,
    thread_counts: list[int],
    tmp_dir: str | Path,
    device: str,
) -> dict[int, float]:
    """Run training for each thread count; return {n_threads: wall_seconds}."""
    tmp_dir = Path(tmp_dir)
    results: dict[int, float] = {}

    n_epochs = _BENCH_CONFIG["num_epochs"]
    cpu_count = os.cpu_count() or "?"
    print(f"Server: {cpu_count} logical CPUs  |  device: {device}")
    print(f"Thread counts to test: {thread_counts}")
    print(f"Epochs per run: {n_epochs}  (dim=64, batch=1024)")
    print()

    _silence_training_logs()

    for n in thread_counts:
        torch.set_num_threads(n)
        print(f"  threads={n:4d} ...", end="", flush=True)
        t0 = time.perf_counter()
        try:
            run_training(
                data_dir=data_dir,
                out_dir=tmp_dir / f"bench_{n}",
                config=_BENCH_CONFIG.copy(),
                device=device,
                cache_dir=cache_dir,
            )
            elapsed = time.perf_counter() - t0
            results[n] = elapsed
            print(f"  {elapsed:6.1f}s")
        except Exception as exc:
            print(f"  FAILED — {exc}")

    return results


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_results(results: dict[int, float], out_path: str | Path) -> None:
    if not results:
        print("No results to plot.")
        return

    threads = sorted(results)
    times   = [results[n] for n in threads]
    best_n  = min(results, key=results.get)
    baseline = results[threads[0]]
    speedups = [baseline / results[n] for n in threads]
    # Ideal linear speedup capped at the max observed speedup
    ideal = [min(n / threads[0], max(speedups)) for n in threads]

    fig, (ax_time, ax_speed) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"Thread benchmark — RotatE  ({_BENCH_CONFIG['num_epochs']} epochs, "
        f"dim={_BENCH_CONFIG['embedding_dim']})",
        fontsize=13, fontweight="bold",
    )

    # ── Left panel: wall time ─────────────────────────────────────────────── #
    colors = ["#2ecc71" if n == best_n else "#4a90d9" for n in threads]
    bars = ax_time.bar([str(n) for n in threads], times, color=colors, width=0.6)
    for bar, t in zip(bars, times):
        ax_time.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(times) * 0.01,
            f"{t:.1f}s", ha="center", va="bottom", fontsize=9,
        )
    ax_time.set_xlabel("OMP_NUM_THREADS", fontsize=11)
    ax_time.set_ylabel("Wall time (seconds)", fontsize=11)
    ax_time.set_title("Training time")
    ax_time.set_ylim(0, max(times) * 1.18)
    ax_time.legend(handles=[
        Patch(facecolor="#2ecc71", label=f"Best: {best_n} threads"),
        Patch(facecolor="#4a90d9", label="Other"),
    ], loc="upper right", fontsize=9)
    ax_time.grid(axis="y", alpha=0.25)

    # ── Right panel: speedup ──────────────────────────────────────────────── #
    ax_speed.plot(
        [str(n) for n in threads], speedups,
        "o-", color="#4a90d9", linewidth=2, markersize=7, label="Actual",
    )
    ax_speed.plot(
        [str(n) for n in threads], ideal,
        "--", color="#bbb", linewidth=1.5, label="Ideal (linear)",
    )
    best_idx = threads.index(best_n)
    ax_speed.axvline(
        x=best_idx, color="#2ecc71", linestyle=":", linewidth=2, alpha=0.9,
        label=f"Best: {best_n} threads",
    )
    for i, (n, s) in enumerate(zip(threads, speedups)):
        ax_speed.annotate(
            f"{s:.2f}×", xy=(i, s),
            xytext=(0, 8), textcoords="offset points",
            ha="center", fontsize=8,
        )
    ax_speed.set_xlabel("OMP_NUM_THREADS", fontsize=11)
    ax_speed.set_ylabel(f"Speedup vs {threads[0]} thread", fontsize=11)
    ax_speed.set_title("Speedup")
    ax_speed.set_ylim(bottom=0)
    ax_speed.legend(fontsize=9)
    ax_speed.grid(axis="y", alpha=0.25)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")

    # ── Print recommendation ──────────────────────────────────────────────── #
    print(f"\nResults:")
    for n in threads:
        marker = " ← best" if n == best_n else ""
        print(f"  {n:4d} threads: {results[n]:.1f}s{marker}")

    print(f"\nPlot saved: {out_path}")
    print(f"\nRecommendation: set OMP_NUM_THREADS={best_n}")
    print(f"  OMP_NUM_THREADS={best_n} uv run python -m src.pykeen.train \\")
    print(f"      --data-dir biocypher_out/ --cache-dir .cache/triples ...")


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Benchmark OMP_NUM_THREADS for RotatE training.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--data-dir",  required=True,
                   help="BioCypher output dir (biocypher_out/).")
    p.add_argument("--cache-dir", default=".cache/triples",
                   help="TriplesFactory cache dir — strongly recommended.")
    p.add_argument("--out",       default="benchmark_threads.png",
                   help="Path for the output PNG chart.")
    p.add_argument("--tmp-dir",   default="/tmp/bench_threads",
                   help="Temporary dir for per-run model outputs.")
    p.add_argument("--max-threads", type=int, default=None,
                   help="Max thread count to test (default: os.cpu_count()).")
    p.add_argument("--device",    default="cpu",
                   choices=["cpu", "cuda", "auto"],
                   help="Device for training (use cpu to benchmark threads).")
    return p.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )

    args = _parse_args()
    max_n = args.max_threads or (os.cpu_count() or 1)
    counts = _thread_counts(max_n)

    results = run_benchmark(
        data_dir=args.data_dir,
        cache_dir=args.cache_dir,
        thread_counts=counts,
        tmp_dir=args.tmp_dir,
        device=args.device,
    )

    plot_results(results, out_path=args.out)
