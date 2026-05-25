"""
run_rgcn.py - CLI entry point for the R-GCN layer ablation study.

Trains one R-GCN encoder per ablation condition, then evaluates multi-label
phenotype-class prediction AUPRC with a fixed gene split and shared class list.

Outputs are written under:
    results/ml/{species}/rgcn/{condition}/

Scientific question answered
-----------------------------
Can molecular graph neighbourhood information (GO, expression, PPI) produce
gene embeddings that are predictive of phenotype-class membership?
Which layers are most informative, and are these patterns species-conserved?

Difference from KGE track
--------------------------
The R-GCN track uses message passing to aggregate neighbourhood context into
gene embeddings, then applies a downstream LogisticRegression classifier.
Phenotype edges are NEVER in the encoder graph - they are labels only.
The KGE track scores triples directly via embedding interactions. See run_kge.py.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from src.machine_learning.common.data_loading import resolve_biocypher_dir
from src.machine_learning.rgcn.ablation import (
    _FAST_PRESET,
    HUMAN_CONDITIONS,
    MOUSE_CONDITIONS,
    run_ablation,
)
from src.machine_learning.rgcn.config import (
    DEFAULT_CONFIG,
    load_config_from_yaml,
    target_relation_for_species,
)

_SPECIES_CONDITIONS = {'human': HUMAN_CONDITIONS, 'mouse': MOUSE_CONDITIONS}
_DEFAULT_CONFIGS = {
    'human': Path('configs/machine_learning/rgcn_human.yaml'),
    'mouse': Path('configs/machine_learning/rgcn_mouse.yaml'),
}

# ── Argument type helpers ─────────────────────────────────────────────────────


def _positive_int(value: str) -> int:
    """
    Convert a string to a positive integer; raise ArgumentTypeError if non-positive.
    """
    n = int(value)
    if n <= 0:
        raise argparse.ArgumentTypeError(f'must be a positive integer, got {value}')
    return n


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
    Parse and return CLI arguments for the R-GCN ablation entry point.

    Returns:
        Parsed argparse.Namespace with all CLI arguments.
    """
    all_conditions = sorted(set(HUMAN_CONDITIONS.keys()) | set(MOUSE_CONDITIONS.keys()))
    p = argparse.ArgumentParser(
        description='R-GCN layer ablation study for gene–phenotype prediction.',
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
        'data/biocypher/<species>_* directory.',
    )
    p.add_argument(
        '--out-dir',
        default=None,
        help='Root output dir. Defaults to results/ml/<species>/rgcn/.',
    )
    p.add_argument(
        '--config',
        default=None,
        help='YAML config file (e.g. configs/machine_learning/rgcn_human.yaml). '
        'Defaults to configs/machine_learning/rgcn_<species>.yaml. '
        'CLI flags override YAML values.',
    )
    p.add_argument(
        '--fast',
        action='store_true',
        help=f'Smoke-test preset: epochs={_FAST_PRESET["rgcn_epochs"]}, '
        f'dim={_FAST_PRESET["dim"]}.',
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
        '--rgcn-epochs',
        type=int,
        default=None,
        help='Number of R-GCN training epochs per condition.',
    )
    p.add_argument(
        '--dim', type=int, default=None, help='R-GCN embedding dimensionality.'
    )
    p.add_argument(
        '--seed',
        type=int,
        default=None,
        help='Random seed for gene split, encoder training, and classifier.',
    )
    p.add_argument(
        '--parallel',
        action='store_true',
        help='Run conditions in parallel (spawn-based subprocesses).',
    )
    p.add_argument(
        '--max-workers',
        type=_positive_int,
        default=None,
        help='Max parallel workers. Defaults to number of selected conditions.',
    )
    p.add_argument(
        '--threads-per-worker',
        type=_positive_int,
        default=None,
        help='OMP_NUM_THREADS per worker. Defaults to cpu_count // max_workers.',
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
    Entry point for the R-GCN layer ablation study.

    Parses CLI arguments, resolves configuration, runs all selected ablation
    conditions via run_ablation(), and prints a summary table to stdout.
    """
    args = _parse_args()

    # Load YAML config. If --config is omitted, use the species default YAML;
    # CLI flags still take precedence over YAML values.
    yaml_cfg: dict = {}
    config_path = (
        Path(args.config) if args.config else _DEFAULT_CONFIGS[args.species or 'human']
    )
    try:
        yaml_cfg = load_config_from_yaml(config_path)
    except FileNotFoundError as exc:
        print(f'ERROR: {exc}')
        sys.exit(1)

    species = _cfg_value(args, yaml_cfg, 'species', 'human')
    fast = args.fast or bool(yaml_cfg.get('fast', False))
    rgcn_epochs = _cfg_value(
        args, yaml_cfg, 'rgcn_epochs', DEFAULT_CONFIG['rgcn_epochs']
    )
    dim = _cfg_value(args, yaml_cfg, 'dim', DEFAULT_CONFIG['dim'])
    if fast:
        rgcn_epochs = _FAST_PRESET['rgcn_epochs']
        dim = _FAST_PRESET['dim']

    seed = _cfg_value(args, yaml_cfg, 'seed', DEFAULT_CONFIG['seed'])
    train_frac = yaml_cfg.get('train_frac', DEFAULT_CONFIG['train_frac'])
    lr = yaml_cfg.get('lr', DEFAULT_CONFIG['lr'])
    batch_size = yaml_cfg.get('batch_size', DEFAULT_CONFIG['batch_size'])
    num_negs_per_pos = yaml_cfg.get(
        'num_negs_per_pos',
        DEFAULT_CONFIG['num_negs_per_pos'],
    )
    classifier_max_iter = yaml_cfg.get(
        'classifier_max_iter',
        DEFAULT_CONFIG['classifier_max_iter'],
    )
    classifier_C = yaml_cfg.get('classifier_C', DEFAULT_CONFIG['classifier_C'])
    parallel = args.parallel or bool(yaml_cfg.get('parallel', False))
    max_workers = _cfg_value(args, yaml_cfg, 'max_workers')
    threads_per_worker = _cfg_value(args, yaml_cfg, 'threads_per_worker')
    conditions_arg = args.conditions or yaml_cfg.get('conditions')
    target_relation = target_relation_for_species(
        species,
        _cfg_value(args, yaml_cfg, 'target_relation'),
    )

    data_dir_arg = _cfg_value(args, yaml_cfg, 'data_dir')
    data_dir = resolve_biocypher_dir(data_dir_arg, species=species)
    _check_data_dir(data_dir, species)

    out_dir = (
        Path(_cfg_value(args, yaml_cfg, 'out_dir'))
        if _cfg_value(args, yaml_cfg, 'out_dir')
        else Path(f'results/ml/{species}/rgcn')
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    selected = _select_conditions(species, conditions_arg)

    print('=' * 70)
    print('R-GCN LAYER ABLATION STUDY')
    print('=' * 70)
    print(f'Config: {config_path.resolve()}')
    print(f'Species: {species}')
    print(f'Target relation: {target_relation}')
    print(f'Data dir: {data_dir.resolve()}')
    print(f'Out dir: {out_dir.resolve()}')
    print(f'Epochs: {rgcn_epochs}  |  Dim: {dim}  |  Seed: {seed}')
    print(f'LR: {lr}  |  Batch: {batch_size}  |  Negs/pos: {num_negs_per_pos}')
    print(f'Fast mode: {fast}')
    print(f'Parallel: {parallel}')
    print(f'Conditions: {", ".join(selected.keys())}')
    print('-' * 70)
    print('NOTE: Phenotype edges are labels only - never in the encoder graph.')
    print('      Gene embeddings are trained via link reconstruction on context edges.')
    print('-' * 70)

    df = run_ablation(
        data_dir=str(data_dir),
        out_dir=str(out_dir),
        conditions=selected,
        species=species,
        target_relation=target_relation,
        rgcn_epochs=rgcn_epochs,
        dim=dim,
        seed=seed,
        train_frac=train_frac,
        lr=lr,
        batch_size=batch_size,
        num_negs_per_pos=num_negs_per_pos,
        classifier_max_iter=classifier_max_iter,
        classifier_C=classifier_C,
        parallel=parallel,
        max_workers=max_workers,
        threads_per_worker=threads_per_worker,
    )

    print()
    print('=' * 70)
    print('R-GCN ABLATION SUMMARY')
    print('=' * 70)
    cols = [
        'condition',
        'mean_auprc',
        'mean_baseline',
        'improvement_x',
        'n_entities',
        'n_context_edges',
        'elapsed_s',
    ]
    show = [c for c in cols if c in df.columns]
    print(df[show].to_string(index=False, float_format='%.4f'))
    print()
    print(f'Full results → {(out_dir / "comparison_summary.csv").resolve()}')
    print(f'Per-condition logs → {out_dir}/<condition>/worker.log')
    print('=' * 70)


if __name__ == '__main__':
    main()
