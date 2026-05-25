"""
src/machine_learning/rgcn/config.py - R-GCN experiment configuration.

Provides Python defaults and YAML-based config loading.
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# ── Target relation names ─────────────────────────────────────────────────────

PHENOTYPE_REL_HUMAN = 'HumanGeneHasMpTopTerm'
PHENOTYPE_REL_MOUSE = 'MouseGeneHasMpTopTerm'

SPECIES_TARGET_RELATIONS = {
    'human': PHENOTYPE_REL_HUMAN,
    'mouse': PHENOTYPE_REL_MOUSE,
}


def target_relation_for_species(
    species: str = 'human',
    target_relation: str | None = None,
) -> str:
    """
    Resolve the R-GCN target gene→phenotype relation used as labels.

    Args:
        species: Species name used when *target_relation* is not explicitly set.
        target_relation: Explicit relation name. When supplied, it wins.

    Returns:
        Relation :TYPE string for gene→MP top-term prediction labels.
    """
    if target_relation:
        return target_relation

    species_key = species.lower()
    try:
        return SPECIES_TARGET_RELATIONS[species_key]
    except KeyError as exc:
        raise ValueError(
            f'Unknown species {species!r}. '
            f'Choose from: {sorted(SPECIES_TARGET_RELATIONS)}'
        ) from exc


# ── Default config ────────────────────────────────────────────────────────────

DEFAULT_CONFIG: dict = {
    'dim': 128,
    'rgcn_epochs': 100,
    'lr': 1e-3,
    'batch_size': 4096,
    'num_negs_per_pos': 5,
    'train_frac': 0.80,
    'seed': 42,
    'classifier_max_iter': 1000,
    'classifier_C': 1.0,
}


def load_config_from_yaml(yaml_path: str | Path) -> dict:
    """
    Load an R-GCN experiment config from a YAML file, merged over Python defaults.

    The YAML may contain any subset of DEFAULT_CONFIG keys plus additional
    fields such as ``species``, ``conditions``, ``data_dir``, and ``out_dir``.
    Unknown keys pass through unmodified to allow experiment-specific extensions.

    Args:
        yaml_path: Path to the YAML config file.

    Returns:
        Dict with YAML values merged over DEFAULT_CONFIG (YAML wins on conflict).

    Raises:
        FileNotFoundError: If the YAML file does not exist.
    """
    import yaml

    yaml_path = Path(yaml_path)
    if not yaml_path.exists():
        raise FileNotFoundError(f'R-GCN config file not found: {yaml_path}')

    with yaml_path.open() as fh:
        overrides = yaml.safe_load(fh) or {}

    merged = {**DEFAULT_CONFIG, **overrides}
    logger.info('Loaded R-GCN config from %s', yaml_path)
    return merged
