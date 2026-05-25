"""
src/machine_learning/kge/config.py - PyKEEN KGE experiment configuration.

Provides Python defaults and YAML-based config loading.
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# ── Target relation names ─────────────────────────────────────────────────────

# Exact :TYPE strings produced by BioCypher after quote-stripping.
GENE_PHENOTYPE_RELATION_HUMAN = 'HumanGeneHasMpTopTerm'
GENE_PHENOTYPE_RELATION_MOUSE = 'MouseGeneHasMpTopTerm'

# Backward-compatible alias - human was the original species.
GENE_PHENOTYPE_RELATION = GENE_PHENOTYPE_RELATION_HUMAN

SPECIES_TARGET_RELATIONS = {
    'human': GENE_PHENOTYPE_RELATION_HUMAN,
    'mouse': GENE_PHENOTYPE_RELATION_MOUSE,
}


def target_relation_for_species(
    species: str = 'human',
    target_relation: str | None = None,
) -> str:
    """
    Resolve the KGE target gene→phenotype relation.

    Args:
        species: Species name used when *target_relation* is not explicitly set.
        target_relation: Explicit relation name. When supplied, it wins.

    Returns:
        Relation :TYPE string for gene→MP top-term prediction.
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


CSV_SEP = '\t'

# ── Default config ────────────────────────────────────────────────────────────

DEFAULT_CONFIG: dict = {
    'model': 'RotatE',
    'embedding_dim': 128,
    'batch_size': 2048,
    'lr': 1e-3,
    'num_negs_per_pos': 64,
    'train_frac': 0.80,
    'val_frac': 0.10,
    'random_seed': 42,
    'eval_batch_size': 128,
    'checkpoint_frequency': 25,
    'target_steps': 80_000,
    'es_patience': 10,
    'es_frequency': 10,
    'es_delta': 2e-3,
}


def load_config_from_yaml(yaml_path: str | Path) -> dict:
    """
    Load a KGE experiment config from a YAML file, merged over Python defaults.

    The YAML may contain any subset of DEFAULT_CONFIG keys plus additional
    fields such as ``species`` and ``conditions``. Unknown keys pass through
    unmodified to allow experiment-specific extensions.

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
        raise FileNotFoundError(f'KGE config file not found: {yaml_path}')

    with yaml_path.open() as fh:
        overrides = yaml.safe_load(fh) or {}

    merged = {**DEFAULT_CONFIG, **overrides}
    logger.info('Loaded KGE config from %s', yaml_path)
    return merged
