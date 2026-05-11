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

    # Quick smoke-test with a tiny config (1 epoch, dim=16):
    python benchmark_threads.py \\
        --data-dir biocypher_out/human/ \\
        --threads 8 16 32 \\
        --epochs 1 --dim 16 --batch 512 --negs 8 \\
        --timeout 300 \\
        --out benchmark_threads.png

    # Preview commands without running anything:
    python benchmark_threads.py \\
        --data-dir biocypher_out/human/ \\
        --threads 8 16 32 \\
        --dry-run
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# When invoked as `python src/pykeen/benchmark_threads.py`, Python puts
# `src/pykeen` on sys.path[0], not the project root, so `src.pykeen.train`
# is unfindable.  Compute the project root from __file__ and prepend it.
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

logger = logging.getLogger(__name__)

# ── Benchmark training config ─────────────────────────────────────────────────
# These defaults produce runs long enough to measure throughput while completing
# in a few minutes.  Pass CLI overrides (--epochs, --dim, --batch, --negs) to
# use a smaller smoke-test config; see DESIGN NOTES at the bottom of this file.
_BENCH_EPOCHS = 5
_BENCH_DIM = 64
_BENCH_BATCH = 1024
_BENCH_NEGS = 32

# ── Threading env var names ───────────────────────────────────────────────────
# All of these cap the CPU thread pool at different layers.  Setting them
# consistently avoids the oversubscription that silently kills throughput.
_THREAD_ENV_VARS = [
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",  # Accelerate (macOS)
    "TF_NUM_INTRAOP_THREADS",  # TensorFlow, no-op for PyKEEN but harmless
    "PYTORCH_NUM_THREADS",
]


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


# ── Preflight checks ──────────────────────────────────────────────────────────


def preflight_checks(
    data_dir: Path,
    cache_dir: Path | None,
    train_module: str,
) -> bool:
    """Run sanity checks before any subprocess is launched.

    Returns True if all checks pass, False if any blocker is found.
    Prints a summary even on success so the log shows the data layout.
    """
    ok = True

    print("\n" + "=" * 68)
    print("PRE-FLIGHT CHECKS")
    print("=" * 68)

    # ── 1. data_dir ───────────────────────────────────────────────────────── #
    if not data_dir.exists():
        print(f"  [FAIL] data_dir does not exist: {data_dir}")
        ok = False
    elif not data_dir.is_dir():
        print(f"  [FAIL] data_dir is a file, expected a directory: {data_dir}")
        ok = False
    else:
        print(f"  [ OK ] data_dir exists: {data_dir.resolve()}")
        csv_files = sorted(data_dir.glob("*.csv"))
        part_files = [f for f in csv_files if "-part" in f.name]
        header_files = [f for f in csv_files if "-header" in f.name]
        print(f"         CSV files total : {len(csv_files)}")
        print(f"         header files    : {len(header_files)}")
        print(f"         data-part files : {len(part_files)}")
        if not csv_files:
            print("  [WARN] No CSV files found — is the BioCypher pipeline done?")
            ok = False
        else:
            # Show approx total data size so we know if CSV parsing will be slow
            total_bytes = sum(f.stat().st_size for f in csv_files)
            print(f"         total CSV size  : {total_bytes / 1_048_576:.1f} MB")
            if total_bytes > 500 * 1_048_576:
                print(
                    "  [WARN] Large dataset (>500 MB) — CSV parsing will be slow "
                    "without a warm --cache-dir.  First run per thread count will "
                    "take extra time; subsequent runs should be much faster."
                )

    # ── 2. cache_dir ──────────────────────────────────────────────────────── #
    if cache_dir is not None:
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)
            print(f"  [ OK ] cache_dir ready  : {cache_dir.resolve()}")
            cached = list(cache_dir.iterdir())
            if cached:
                print(f"         cached files    : {len(cached)} (TriplesFactory cache present)")
            else:
                print(
                    "         cached files    : 0  "
                    "(first run will parse CSVs and write cache)"
                )
        except OSError as exc:
            print(f"  [FAIL] cannot create cache_dir {cache_dir}: {exc}")
            ok = False
    else:
        print(
            "  [WARN] --cache-dir not set.  Every subprocess will re-parse all "
            "CSV files from scratch.  This can add 1–5 minutes per run.  "
            "Pass --cache-dir .cache/triples to enable caching."
        )

    # ── 3. train module importability ─────────────────────────────────────── #
    # We check at the *spec* level (find the module without importing it) to
    # avoid loading PyTorch/PyKEEN here, which would defeat the purpose of the
    # subprocess isolation pattern.
    #
    # find_spec() raises ModuleNotFoundError (not just returns None) when a
    # parent package exists on disk but is not importable.  Catch it so the
    # preflight summary is never an unhandled traceback.
    try:
        spec = importlib.util.find_spec(train_module)
        _find_spec_exc: BaseException | None = None
    except ModuleNotFoundError as _exc:
        spec = None
        _find_spec_exc = _exc

    if spec is not None:
        print(f"  [ OK ] module found     : {train_module}  →  {spec.origin}")
    else:
        # Try to give a more helpful diagnosis
        parts = train_module.split(".")
        parent = ".".join(parts[:-1])
        try:
            parent_spec = importlib.util.find_spec(parent) if parent else None
        except ModuleNotFoundError:
            parent_spec = None
        exc_detail = f" ({_find_spec_exc})" if _find_spec_exc else ""
        if parent_spec is None:
            print(
                f"  [WARN] module not found : {train_module}{exc_detail}\n"
                f"         Parent package '{parent}' is also missing from "
                f"sys.path.  Make sure you run this script from the project "
                f"root and that the package is installed or on PYTHONPATH.\n"
                f"         sys.path:\n"
                + "\n".join(f"           {p}" for p in sys.path)
            )
        else:
            print(
                f"  [WARN] module not found : {train_module}{exc_detail}\n"
                f"         Parent package '{parent}' exists at {parent_spec.origin}, "
                f"but the 'train' sub-module was not found.  "
                f"Check the module path."
            )
        # Not a hard blocker — maybe uv resolves it differently in subprocesses
        # but warn the user clearly.

    # ── 4. Python executable ──────────────────────────────────────────────── #
    print(f"  [ OK ] Python           : {sys.executable}  ({sys.version})")
    print(f"  [ OK ] cwd              : {Path.cwd()}")

    print("=" * 68 + "\n")
    return ok


# ── Per-run logging setup ──────────────────────────────────────────────────────


def _make_run_logger(log_path: Path, n_threads: int) -> logging.Logger:
    """Return a dedicated Logger that writes to *log_path*."""
    run_logger = logging.getLogger(f"bench.threads.{n_threads}")
    run_logger.setLevel(logging.DEBUG)
    run_logger.propagate = False  # do not bubble to root logger

    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(
        logging.Formatter(
            "%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%H:%M:%S"
        )
    )
    run_logger.addHandler(fh)
    return run_logger


# ── Build subprocess command + env ────────────────────────────────────────────


def _build_cmd_and_env(
    n_threads: int,
    data_dir: Path,
    out_dir: Path,
    cache_dir: Path | None,
    epochs: int,
    dim: int,
    batch: int,
    negs: int,
    train_module: str,
) -> tuple[list[str], dict[str, str]]:
    """Return (cmd, env) for the subprocess call."""
    env = os.environ.copy()

    # Ensure the project root reaches the subprocess so `python -m src.pykeen.train`
    # resolves even when the subprocess inherits a stripped or absent PYTHONPATH.
    _existing_pp = env.get("PYTHONPATH", "")
    _root_str = str(_PROJECT_ROOT)
    if _root_str not in _existing_pp.split(os.pathsep):
        env["PYTHONPATH"] = _root_str + (os.pathsep + _existing_pp if _existing_pp else "")

    for var in _THREAD_ENV_VARS:
        env[var] = str(n_threads)

    cmd = [
        sys.executable, "-m", train_module,
        "--data-dir", str(data_dir.resolve()),
        "--out-dir",  str(out_dir.resolve()),
        "--epochs",   str(epochs),
        "--dim",      str(dim),
        "--batch",    str(batch),
        "--negs",     str(negs),
        "--device",   "cpu",
    ]
    if cache_dir is not None:
        # Always pass an *absolute* cache path so every subprocess uses the same
        # location regardless of their working directory.
        cmd += ["--cache-dir", str(cache_dir.resolve())]

    return cmd, env


# ── Single subprocess run (streaming) ────────────────────────────────────────


def _run_one(
    n_threads: int,
    data_dir: Path,
    out_dir: Path,
    cache_dir: Path | None,
    log_dir: Path,
    timeout: int,
    epochs: int,
    dim: int,
    batch: int,
    negs: int,
    train_module: str,
    dry_run: bool = False,
) -> float | None:
    """Launch one training run in a subprocess with OMP_NUM_THREADS=n_threads.

    Streams stdout/stderr line-by-line to both the console and a dedicated log
    file at  <log_dir>/threads_<n>.log .

    Returns wall-clock seconds, or None on failure.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    log_path = log_dir / f"threads_{n_threads}.log"
    run_log = _make_run_logger(log_path, n_threads)

    cmd, env = _build_cmd_and_env(
        n_threads, data_dir, out_dir, cache_dir,
        epochs, dim, batch, negs, train_module,
    )

    # ── Header block written to the per-run log ──────────────────────────── #
    run_log.info("=" * 60)
    run_log.info("BENCHMARK RUN — OMP_NUM_THREADS=%d", n_threads)
    run_log.info("=" * 60)
    run_log.info("Command   : %s", " ".join(cmd))
    run_log.info("cwd       : %s", Path.cwd())
    run_log.info("Python    : %s", sys.executable)
    run_log.info("data_dir  : %s", data_dir.resolve())
    run_log.info("out_dir   : %s", out_dir.resolve())
    run_log.info("cache_dir : %s", cache_dir.resolve() if cache_dir else "disabled")
    run_log.info("timeout   : %ds", timeout)
    run_log.info("config    : epochs=%d  dim=%d  batch=%d  negs=%d", epochs, dim, batch, negs)
    run_log.info("--- threading env vars ---")
    for var in _THREAD_ENV_VARS + ["PYTORCH_NUM_THREADS"]:
        run_log.info("  %s = %s", var, env.get(var, "(not set)"))
    run_log.info("--- subprocess stdout/stderr below ---")

    if dry_run:
        print(f"    [dry-run] would execute:\n      {' '.join(cmd)}")
        print(f"    [dry-run] log would be written to: {log_path}")
        run_log.info("[dry-run] skipping execution")
        return None

    t0 = time.perf_counter()
    proc: subprocess.Popen | None = None

    try:
        proc = subprocess.Popen(
            cmd,
            env=env,
            cwd=Path.cwd(),           # inherit parent cwd explicitly
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # merge stderr into stdout
            text=True,
            bufsize=1,                 # line-buffered
        )

        # Stream output line-by-line.  This is the critical change: instead of
        # waiting for the process to exit (or timeout), we see every log line
        # the moment the child writes it.  This makes it possible to identify
        # exactly where the child process is hanging.
        assert proc.stdout is not None
        deadline = t0 + timeout

        for line in proc.stdout:
            line = line.rstrip("\n")
            run_log.info("%s", line)

            # Also mirror to the parent console so you see progress live.
            # Prefix with thread count so interleaved output is readable.
            print(f"  [t={n_threads}] {line}", flush=True)

            if time.perf_counter() > deadline:
                # We have been reading output past the deadline — the process
                # is still alive.  Break out and terminate it below.
                run_log.warning(
                    "Deadline reached after %.0fs — terminating subprocess",
                    time.perf_counter() - t0,
                )
                break

        # Give the process a moment to exit on its own after the loop ends
        # (either stdout closed naturally, or we broke out on deadline).
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            pass

    except Exception as exc:
        run_log.error("Unexpected error launching subprocess: %s", exc, exc_info=True)
        logger.error("threads=%d: subprocess launch error: %s", n_threads, exc)
        if proc is not None:
            _terminate_proc(proc, n_threads, run_log)
        return None

    # ── Check whether deadline was exceeded ──────────────────────────────── #
    elapsed = time.perf_counter() - t0

    if proc.poll() is None:
        # Process is still alive — timeout scenario
        run_log.error(
            "TIMEOUT after %.1fs (limit=%ds) — sending SIGTERM, then SIGKILL",
            elapsed, timeout,
        )
        logger.error(
            "Run with threads=%d timed out after %.0fs", n_threads, elapsed
        )
        _terminate_proc(proc, n_threads, run_log)
        return None

    rc = proc.returncode
    run_log.info("Process exited with returncode=%d  elapsed=%.1fs", rc, elapsed)

    if rc != 0:
        run_log.error(
            "Non-zero returncode=%d — check log above for the first ERROR line",
            rc,
        )
        logger.error(
            "Run with threads=%d failed (rc=%d) — see %s", n_threads, rc, log_path
        )
        return None

    run_log.info("SUCCESS  elapsed=%.1fs", elapsed)
    return elapsed


def _terminate_proc(
    proc: subprocess.Popen,
    n_threads: int,
    run_log: logging.Logger,
) -> None:
    """Gracefully terminate *proc*: SIGTERM → wait 15s → SIGKILL."""
    if proc.poll() is not None:
        return  # already dead
    run_log.warning("Sending SIGTERM to pid=%d", proc.pid)
    proc.terminate()
    try:
        proc.wait(timeout=15)
        run_log.info("Process exited after SIGTERM")
    except subprocess.TimeoutExpired:
        run_log.warning("Process did not exit after SIGTERM; sending SIGKILL")
        proc.kill()
        proc.wait()
        run_log.info("Process killed")


# ── Benchmark loop ────────────────────────────────────────────────────────────


def run_benchmark(
    data_dir: Path,
    thread_counts: list[int],
    tmp_root: Path,
    log_dir: Path,
    cache_dir: Path | None = None,
    timeout: int = 3600,
    epochs: int = _BENCH_EPOCHS,
    dim: int = _BENCH_DIM,
    batch: int = _BENCH_BATCH,
    negs: int = _BENCH_NEGS,
    train_module: str = "src.pykeen.train",
    dry_run: bool = False,
) -> dict[int, float]:
    """Run one subprocess per thread count; return {n_threads: wall_seconds}."""
    results: dict[int, float] = {}

    print(f"CPU count   : {os.cpu_count()}")
    print(f"Threads     : {thread_counts}")
    print(f"Config      : epochs={epochs}  dim={dim}  batch={batch}  negs={negs}")
    print(f"Timeout     : {timeout}s per run")
    print(f"Cache dir   : {cache_dir.resolve() if cache_dir else 'disabled (slow)'}")
    print(f"Log dir     : {log_dir.resolve()}")
    print(f"Dry run     : {dry_run}")
    print()

    # Oversubscription warning — on large servers this is a real concern.
    cpu_count = os.cpu_count() or 1
    for n in thread_counts:
        if n > cpu_count:
            print(
                f"  [WARN] threads={n} > cpu_count={cpu_count}: "
                f"this will oversubscribe the CPU and is likely slower."
            )

    print()

    for n in thread_counts:
        out_dir = tmp_root / f"bench_{n}"
        print(f"  threads={n:4d} ... ", end="", flush=True)
        t_start = time.perf_counter()

        elapsed = _run_one(
            n_threads=n,
            data_dir=data_dir,
            out_dir=out_dir,
            cache_dir=cache_dir,
            log_dir=log_dir,
            timeout=timeout,
            epochs=epochs,
            dim=dim,
            batch=batch,
            negs=negs,
            train_module=train_module,
            dry_run=dry_run,
        )

        if elapsed is None:
            if not dry_run:
                print(f"  FAILED  (see {log_dir}/threads_{n}.log)")
        else:
            results[n] = elapsed
            print(f"  {elapsed:6.1f}s")

    return results


# ── Plotting ──────────────────────────────────────────────────────────────────


def plot_results(
    results: dict[int, float],
    out_path: Path,
    epochs: int,
    dim: int,
    batch: int,
) -> None:
    if not results:
        print("No results to plot.")
        return

    threads = sorted(results)
    times = [results[n] for n in threads]
    best_n = min(results, key=results.__getitem__)
    baseline = results[threads[0]]
    speedups = [baseline / results[n] for n in threads]
    ideal = [min(n / threads[0], max(speedups)) for n in threads]

    fig, (ax_t, ax_s) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"OMP_NUM_THREADS benchmark — RotatE  "
        f"({epochs} epochs, dim={dim}, batch={batch})",
        fontsize=13, fontweight="bold",
    )

    # ── Left: wall time ────────────────────────────────────────────────────── #
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
    ax_t.legend(
        handles=[
            Patch(facecolor="#2ecc71", label=f"Best: {best_n} threads"),
            Patch(facecolor="#4a90d9", label="Other"),
        ],
        fontsize=9,
    )
    ax_t.grid(axis="y", alpha=0.25)

    # ── Right: speedup ─────────────────────────────────────────────────────── #
    xs = [str(n) for n in threads]
    ax_s.plot(xs, speedups, "o-", color="#4a90d9", lw=2, ms=7, label="Actual")
    ax_s.plot(xs, ideal, "--", color="#bbb", lw=1.5, label="Ideal (linear)")
    ax_s.axvline(
        x=str(best_n), color="#2ecc71",
        linestyle=":", lw=2, alpha=0.9, label=f"Best: {best_n} threads",
    )
    for i, (n, s) in enumerate(zip(threads, speedups)):
        ax_s.annotate(
            f"{s:.2f}×",
            xy=(i, s),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            fontsize=8,
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

    # ── Console summary ─────────────────────────────────────────────────────── #
    print("\nResults:")
    for n in threads:
        marker = "  ← best" if n == best_n else ""
        print(f"  {n:4d} threads:  {results[n]:.1f}s{marker}")

    print(f"\nPlot saved: {out_path}")
    print(f"\nRecommendation:")
    print(f"  OMP_NUM_THREADS={best_n} MKL_NUM_THREADS={best_n} \\")
    print(f"  uv run python -m src.pykeen.train --data-dir biocypher_out/human/ ...")

    json_path = out_path.with_suffix(".json")
    json_path.write_text(
        json.dumps(
            {"thread_counts": threads, "wall_seconds": times, "best_n_threads": best_n},
            indent=2,
        )
    )
    print(f"Raw data  : {json_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Benchmark OMP_NUM_THREADS for RotatE (subprocess-based).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--data-dir", required=True,
        help="BioCypher output dir, e.g. biocypher_out/human/"
    )
    p.add_argument(
        "--out", default="benchmark_threads.png",
        help="Output PNG path."
    )
    p.add_argument(
        "--threads", type=int, nargs="+", default=None,
        help="Explicit thread counts to test.  "
             "Default: powers of 2 up to os.cpu_count().",
    )

    # ── Per-run training config ──────────────────────────────────────────── #
    p.add_argument("--epochs", type=int, default=_BENCH_EPOCHS,
                   help="Training epochs per subprocess run.")
    p.add_argument("--dim",    type=int, default=_BENCH_DIM,
                   help="Embedding dimension (must be even for RotatE).")
    p.add_argument("--batch",  type=int, default=_BENCH_BATCH,
                   help="Training batch size.")
    p.add_argument("--negs",   type=int, default=_BENCH_NEGS,
                   help="Negative samples per positive triple.")

    # ── Infrastructure ───────────────────────────────────────────────────── #
    p.add_argument(
        "--cache-dir", default=".cache/triples",
        help="TriplesFactory cache dir.  Strongly recommended: without it "
             "each subprocess re-parses all CSV files (~1–5 min each).",
    )
    p.add_argument(
        "--tmp-dir", default=None,
        help="Temp dir for per-run outputs (auto-created if omitted).",
    )
    p.add_argument(
        "--log-dir", default="benchmark_logs",
        help="Directory for per-thread-count log files.",
    )
    p.add_argument(
        "--timeout", type=int, default=3600,
        help="Per-subprocess wall-clock timeout in seconds.",
    )
    p.add_argument(
        "--train-module", default="src.pykeen.train",
        help="Python module path for the training entry point.",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="Print commands and environment without executing anything.",
    )
    p.add_argument(
        "--skip-preflight", action="store_true",
        help="Skip pre-flight sanity checks (data_dir, cache, module import).",
    )
    return p.parse_args()


# ═══════════════════════════════════════════════════════════════════════════════
# DESIGN NOTES
# ═══════════════════════════════════════════════════════════════════════════════
#
# WHY EVERY RUN TIMED OUT (original script)
# ──────────────────────────────────────────
# The original script used subprocess.run(..., stdout=PIPE, stderr=STDOUT,
# timeout=3600).  This call blocks silently until the child exits or the 60-min
# timer fires; you see no output until then.  The result: five silent 60-minute
# runs, no diagnostic information, 5 × FAILED.
#
# Possible actual bottlenecks (in rough probability order):
#
#   1. CSV parsing / TriplesFactory build on first run.
#      With a large BioCypher graph (> 500 MB of CSVs) the pandas read loop in
#      loader.py can take 5–20 minutes.  With stdout buffered you never see
#      "Loading TriplesFactory from cache" or "Found N edge types", so you
#      cannot tell whether you are stuck here.
#      FIX: use --cache-dir and make sure it points to an absolute path so every
#      subprocess writes/reads from the same location; the new script resolves
#      cache_dir to an absolute path before passing it to the subprocess.
#
#   2. Thread oversubscription with high thread counts.
#      On a 255-core server, PyTorch's OpenMP thread pool can spawn hundreds of
#      OS threads.  When OMP_NUM_THREADS=128 and you have a large batch, BLAS
#      and OpenMP can create 128 × N worker threads.  The OS scheduler thrashes
#      and the process appears frozen.
#      FIX: the new script sets OMP, MKL, OPENBLAS, NUMEXPR, and
#      PYTORCH_NUM_THREADS all to the same value.
#
#   3. Relative cache_dir path.
#      The original script passed --cache-dir .cache/triples which is a
#      *relative* path.  If the subprocess's cwd differs from the parent, each
#      run creates a fresh cache in a different location and re-parses the CSVs
#      every single time.  The new _build_cmd_and_env() always resolves to an
#      absolute path via cache_dir.resolve().
#
#   4. Module path or PYTHONPATH issue.
#      If `src.pykeen.train` is not importable in the subprocess (e.g. the
#      package is not installed and the cwd is different), the subprocess exits
#      immediately with a non-zero code — but with stdout buffered you never see
#      the ModuleNotFoundError.  The new script streams output so you see this
#      immediately.
#
# SMOKE-TEST BEFORE A FULL BENCHMARK
# ────────────────────────────────────
# Before committing to 5 epochs × 5 thread counts (potentially 5+ hours), run
# a single sanity check with the smallest possible config:
#
#   python benchmark_threads.py \
#       --data-dir biocypher_out/ \
#       --threads 8 \
#       --epochs 1 --dim 16 --batch 256 --negs 4 \
#       --timeout 600 \
#       --out /tmp/smoke.png
#
# This should complete in under 10 minutes.  If it times out, check the
# benchmark_logs/threads_8.log file — the streaming output will show exactly
# which step (CSV read, cache build, model init, epoch 1) is the bottleneck.
#
# THREAD COUNT SELECTION FOR A 255-CPU SERVER
# ────────────────────────────────────────────
# RotatE on CPU is primarily memory-bandwidth bound (embedding look-ups +
# complex multiplications), not compute bound.  On modern multi-socket NUMA
# servers the sweet spot is typically 1–2× the number of physical cores per
# socket, not the total logical core count.  Suggested starting points:
#   --threads 4 8 16 32 64
# Testing 128+ on a 255-core server is valid but expect diminishing or negative
# returns past the per-socket core count due to NUMA-remote memory access.
#
# ═══════════════════════════════════════════════════════════════════════════════


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )

    args = _parse_args()
    counts = args.threads or _default_thread_counts()
    data_dir = Path(args.data_dir)
    cache_dir = Path(args.cache_dir) if args.cache_dir else None
    log_dir = Path(args.log_dir)

    if not args.skip_preflight:
        preflight_ok = preflight_checks(data_dir, cache_dir, args.train_module)
        if not preflight_ok and not args.dry_run:
            print(
                "\nPre-flight checks found problems.  "
                "Fix them before benchmarking, or re-run with --skip-preflight "
                "to proceed anyway."
            )
            sys.exit(1)

    own_tmp = args.tmp_dir is None
    tmp_root = (
        Path(args.tmp_dir) if args.tmp_dir
        else Path(tempfile.mkdtemp(prefix="bench_threads_"))
    )

    try:
        results = run_benchmark(
            data_dir=data_dir,
            thread_counts=counts,
            tmp_root=tmp_root,
            log_dir=log_dir,
            cache_dir=cache_dir,
            timeout=args.timeout,
            epochs=args.epochs,
            dim=args.dim,
            batch=args.batch,
            negs=args.negs,
            train_module=args.train_module,
            dry_run=args.dry_run,
        )
        if not args.dry_run:
            plot_results(
                results,
                out_path=Path(args.out),
                epochs=args.epochs,
                dim=args.dim,
                batch=args.batch,
            )
    finally:
        if own_tmp and tmp_root.exists():
            shutil.rmtree(tmp_root, ignore_errors=True)