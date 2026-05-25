"""
run_kge.py - CLI entry point for the KGE (PyKEEN) layer ablation study.

Trains one KGE model per ablation condition, evaluates AUPRC against a fixed
test set (no contamination), and writes per-condition artifacts under:
    results/ml/{species}/kge/{model_lower}/{condition}/

Scientific question answered
-----------------------------
Can the KGE model score missing gene→phenotype triples?
Which molecular layers (GO, expression, PPI) contribute to link prediction
performance, and are these patterns conserved between human and mouse?

Difference from R-GCN track
----------------------------
The KGE track scores triples directly via embedding dot products.
The R-GCN track uses message passing to build gene embeddings, then applies
a downstream classifier. See run_rgcn.py for the R-GCN entry point.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from src.machine_learning.kge.ablation import (
    HUMAN_CONDITIONS,
    MOUSE_CONDITIONS,
    run_ablation,
)
from src.machine_learning.kge.config import (
    load_config_from_yaml,
    target_relation_for_species,
)
from src.machine_learning.kge.loader import resolve_biocypher_dir
from src.utils.logging_config import setup_logging

_SPECIES_CONDITIONS = {'human': HUMAN_CONDITIONS, 'mouse': MOUSE_CONDITIONS}

# ── Argument type helpers ─────────────────────────────────────────────────────


def _positive_int(value: str) -> int:
    """
    Convert a string to a positive integer; raise ArgumentTypeError if non-positive.
    """
    n = int(value)
    if n <= 0:
        raise argparse.ArgumentTypeError(f'must be a positive integer, got {value}')
    return n


# ── Thread environment helpers ────────────────────────────────────────────────

_THREAD_ENV_VARS = [
    'OMP_NUM_THREADS',
    'MKL_NUM_THREADS',
    'OPENBLAS_NUM_THREADS',
    'NUMEXPR_NUM_THREADS',
    'PYTORCH_NUM_THREADS',
]


def _apply_threads(n: int | None) -> None:
    """
    Set all thread-count environment variables to *n*, or warn if OMP_NUM_THREADS is unset.
    """
    if n is not None:
        for var in _THREAD_ENV_VARS:
            os.environ[var] = str(n)
    elif not os.environ.get('OMP_NUM_THREADS'):
        print(
            f'WARNING: --threads not set and OMP_NUM_THREADS is unset. '
            f'PyTorch will use all {os.cpu_count()} CPU cores, which is '
            'usually suboptimal for KGE training on large servers.\n'
        )


# ── Data directory validation ─────────────────────────────────────────────────


def _check_data_dir(data_dir: Path, species: str) -> None:
    """
    Validate that *data_dir* exists and contains BioCypher CSV header files.

    Prints an error message and exits with code 1 if validation fails.

    Args:
        data_dir: Path to the BioCypher output directory.
        species: Species label used in the error hint command.
    """
    if not data_dir.exists():
        print(f'ERROR: data directory does not exist: {data_dir.resolve()}')
        print('\nBuild the knowledge graph first:')
        print(f'  uv run python run_kg.py --species {species}\n')
        sys.exit(1)

    headers = list(data_dir.glob('*-header.csv'))
    if not headers:
        print(f'ERROR: no BioCypher CSV files found in {data_dir.resolve()}')
        print('Re-run the graph build:')
        print(f'  uv run python run_kg.py --species {species}\n')
        sys.exit(1)

    print(f'Data directory OK - {len(headers)} relation type(s) found in {data_dir}\n')


# ── CLI ───────────────────────────────────────────────────────────────────────


def _parse_args() -> argparse.Namespace:
    """
    Parse and return CLI arguments for the KGE ablation entry point.

    Returns:
        Parsed argparse.Namespace with all CLI arguments.
    """
    all_conditions = sorted(set(HUMAN_CONDITIONS.keys()) | set(MOUSE_CONDITIONS.keys()))
    p = argparse.ArgumentParser(
        description='KGE (PyKEEN) layer ablation study for gene–phenotype prediction.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        '--species',
        default=None,
        choices=['human', 'mouse'],
        help='Species to run. Selects the correct graph and relation names.',
    )
    p.add_argument(
        '--target-relation',
        default=None,
        help='Explicit gene→phenotype target relation. Overrides --species.',
    )
    p.add_argument(
        '--data-dir',
        default=None,
        help='BioCypher output directory. Defaults to the newest '
        'data/biocypher/<species>_* directory produced by run_kg.py.',
    )
    p.add_argument(
        '--out-dir',
        default=None,
        help='Root output dir. Defaults to results/ml/<species>/kge/<model_lower>/.',
    )
    p.add_argument(
        '--config',
        default=None,
        help='YAML config file (e.g. configs/machine_learning/kge_human.yaml). '
        'Values override Python defaults. CLI flags override YAML.',
    )
    p.add_argument(
        '--fast',
        action='store_true',
        help='Smoke-test preset: dim=64, 30 epochs per condition.',
    )
    p.add_argument(
        '--conditions',
        nargs='+',
        default=None,
        choices=all_conditions,
        metavar='CONDITION',
        help=f'Subset of conditions to run. Choices: {", ".join(all_conditions)}.',
    )
    p.add_argument(
        '--threads',
        type=_positive_int,
        default=None,
        help='CPU thread count (sets OMP_NUM_THREADS and related vars).',
    )
    p.add_argument(
        '--model',
        default=None,
        choices=['RotatE', 'ComplEx', 'DistMult', 'TransE'],
        help='KGE model to train for every condition.',
    )
    p.add_argument(
        '--seed',
        type=int,
        default=None,
        help='Random seed for train/validation/test split and training.',
    )
    p.add_argument(
        '--target-steps',
        type=int,
        default=None,
        help='Gradient-step budget per condition. Overrides num_epochs so every '
        'condition trains for the same number of parameter updates.',
    )
    return p.parse_args()


def _cfg_value(
    args: argparse.Namespace,
    cfg: dict,
    attr: str,
    default=None,
    *,
    cfg_key: str | None = None,
):
    """
    Return the effective config value with CLI > YAML > default precedence.

    Args:
        args: Parsed argparse namespace.
        cfg: YAML config dict.
        attr: Attribute name on *args* (also used as YAML key unless *cfg_key* is set).
        default: Fallback if neither CLI nor YAML provides a value.
        cfg_key: Override YAML key when it differs from *attr*.

    Returns:
        The first non-None value from: CLI arg, YAML config, default.
    """
    value = getattr(args, attr)
    if value is not None:
        return value
    return cfg.get(cfg_key or attr, default)


def _select_conditions(
    species: str,
    requested: list[str] | None,
) -> dict[str, set[str] | None]:
    """
    Return the condition dict for *species*, filtered to *requested* if given.

    Args:
        species: 'human' or 'mouse'.
        requested: Optional list of condition names to run. None returns all.

    Returns:
        Dict mapping condition name → set of active relation names (or None).
    """
    species_conditions = _SPECIES_CONDITIONS.get(species, HUMAN_CONDITIONS)
    if not requested:
        return species_conditions

    unknown = sorted(set(requested) - set(species_conditions))
    if unknown:
        print(f'ERROR: unknown condition(s) for {species}: {", ".join(unknown)}')
        print(f'Choose from: {", ".join(sorted(species_conditions))}')
        sys.exit(1)

    return {k: species_conditions[k] for k in requested}


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> None:
    """
    Entry point for the KGE layer ablation study.

    Parses CLI arguments, resolves configuration, runs all selected ablation
    conditions via run_ablation(), and prints a summary table to stdout.
    """
    args = _parse_args()

    # Load YAML first so config can supply species, data_dir, out_dir, etc.
    yaml_config: dict = {}
    if args.config:
        try:
            yaml_config = load_config_from_yaml(args.config)
        except FileNotFoundError as exc:
            print(f'ERROR: {exc}')
            sys.exit(1)

    species = _cfg_value(args, yaml_config, 'species', 'human')
    model = _cfg_value(args, yaml_config, 'model', 'RotatE')
    data_dir_arg = _cfg_value(args, yaml_config, 'data_dir')
    out_dir_arg = _cfg_value(args, yaml_config, 'out_dir')
    target_relation_arg = _cfg_value(args, yaml_config, 'target_relation')
    target_steps = _cfg_value(args, yaml_config, 'target_steps')
    seed = _cfg_value(args, yaml_config, 'seed', cfg_key='random_seed')
    fast = args.fast or bool(yaml_config.get('fast', False))
    conditions_arg = args.conditions or yaml_config.get('conditions')
    threads = _cfg_value(args, yaml_config, 'threads')

    # Thread env must be set before torch is imported anywhere.
    _apply_threads(threads)

    data_dir = resolve_biocypher_dir(data_dir_arg, species=species)
    _check_data_dir(data_dir, species)

    # Resolve out_dir: CLI > YAML > default
    out_dir = (
        Path(out_dir_arg)
        if out_dir_arg
        else Path(f'results/ml/{species}/kge/{model.lower()}')
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    setup_logging(out_dir=out_dir, log_filename='ablation.log')

    target_relation = target_relation_for_species(
        species,
        target_relation_arg,
    )

    selected = _select_conditions(species, conditions_arg)

    n_threads = threads or os.environ.get('OMP_NUM_THREADS') or os.cpu_count()

    print('=' * 70)
    print('KGE LAYER ABLATION STUDY')
    print('=' * 70)
    print(f'Species: {species}')
    print(f'Target relation: {target_relation}')
    print(f'Model: {model}')
    print(f'Data dir: {data_dir.resolve()}')
    print(f'Out dir: {out_dir.resolve()}')
    print(f'Fast mode: {fast}')
    print(f'Seed: {seed}')
    print(f'Threads: {n_threads}')
    print(f'Conditions: {", ".join(selected.keys())}')
    print('-' * 70)

    df = run_ablation(
        data_dir=data_dir,
        out_dir=out_dir,
        conditions=selected,
        species=species,
        target_relation=target_relation,
        fast=fast,
        target_steps=target_steps,
        model=model,
        config={
            **yaml_config,
            **({'random_seed': seed} if seed is not None else {}),
        },
    )

    print()
    print('=' * 70)
    print('KGE ABLATION SUMMARY')
    print('=' * 70)
    cols = ['condition', 'mean_auprc', 'mrr', 'hits_at_1', 'hits_at_3', 'hits_at_10']
    print(
        df[[c for c in cols if c in df.columns]].to_string(
            index=False, float_format='%.4f'
        )
    )
    print()
    print(f'Full results → {(out_dir / "comparison_summary.csv").resolve()}')
    print(f'Log → {(out_dir / "ablation.log").resolve()}')
    print('=' * 70)


if __name__ == '__main__':
    main()
