"""
src/machine_learning/kge/loader.py - PyKEEN TriplesFactory builder for the KGE track.

Wraps common BioCypher CSV loading with TriplesFactory construction, optional
binary caching, and a leakage sanity check for the target relation.

Re-exports resolve_biocypher_dir / discover_edge_types / load_triples from
common/data_loading so that callers can import them from one place.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

from pykeen.triples import TriplesFactory
from src.machine_learning.common.data_loading import (  # noqa: F401 - re-exported for callers
    CSV_SEP,
    _read_type,
    discover_edge_types,
    load_triples,
    resolve_biocypher_dir,
)

logger = logging.getLogger(__name__)

__all__ = [
    'resolve_biocypher_dir',
    'discover_edge_types',
    'load_triples',
    'build_triples_factory',
]


# ── TriplesFactory cache helpers ──────────────────────────────────────────────


def _tf_cache_key(
    include_relations: set[str] | None,
    create_inverse_triples: bool,
    source_files: frozenset[str] = frozenset(),
) -> str:
    """
    Return a short deterministic key encoding loader parameters and source file set.
    """
    payload = json.dumps(
        {
            'rels': sorted(include_relations) if include_relations else None,
            'inv': create_inverse_triples,
            'files': sorted(source_files),
        },
        sort_keys=True,
    )
    return hashlib.md5(payload.encode()).hexdigest()[:10]


def _cache_is_fresh(cache_path: Path, source_dir: Path) -> bool:
    """
    Return True if *cache_path* exists and is newer than every CSV in *source_dir*.
    """
    if not cache_path.exists():
        return False
    cache_mtime = cache_path.stat().st_mtime
    return all(f.stat().st_mtime <= cache_mtime for f in source_dir.glob('*.csv'))


# ── Public API ────────────────────────────────────────────────────────────────


def build_triples_factory(
    bc_dir: str | Path | None,
    include_relations: set[str] | None = None,
    create_inverse_triples: bool = False,
    sep: str = CSV_SEP,
    cache_dir: str | Path | None = None,
    target_relation: str | None = None,
) -> TriplesFactory:
    """
    Build a PyKEEN TriplesFactory from a BioCypher output directory.

    Args:
        bc_dir: BioCypher output dir (None → auto-detect newest human build).
        include_relations: Relation type names to include. None loads all.
        create_inverse_triples: Add (tail, relation_inverse, head) for every triple.
            Can improve tail-prediction quality at the cost of doubling triple count.
        sep: CSV field separator.
        cache_dir: Directory for caching the TriplesFactory binary. Cache is
            invalidated when any source CSV changes or when the CSV file set changes.
        target_relation: Optional gene→phenotype relation to check for.

    Returns:
        TriplesFactory ready for training or evaluation.
    """
    bc_dir = resolve_biocypher_dir(bc_dir)

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        source_files = frozenset(f.name for f in bc_dir.glob('*.csv'))
        key = _tf_cache_key(include_relations, create_inverse_triples, source_files)
        cache_path = cache_dir / f'tf_{key}'
        if _cache_is_fresh(cache_path, bc_dir):
            logger.info('Loading TriplesFactory from cache: %s', cache_path)
            tf = TriplesFactory.from_path_binary(cache_path)
            logger.info(
                'TriplesFactory ready (cached) - %d entities | %d relations | %d triples',
                tf.num_entities,
                tf.num_relations,
                tf.num_triples,
            )
            return tf
        logger.info('Cache miss (key=%s) - building from CSV files', key)

    triples = load_triples(bc_dir, include_relations=include_relations, sep=sep)

    tf = TriplesFactory.from_labeled_triples(
        triples=triples,
        create_inverse_triples=create_inverse_triples,
    )

    logger.info(
        'TriplesFactory ready - %d entities | %d relations | %d triples',
        tf.num_entities,
        tf.num_relations,
        tf.num_triples,
    )

    if cache_dir is not None:
        tf.to_path_binary(cache_path)
        logger.info('TriplesFactory cached to: %s', cache_path)

    if target_relation is not None and target_relation not in tf.relation_to_id:
        logger.warning(
            "Primary relation '%s' not found in the loaded data!\n"
            'Available relations: %s',
            target_relation,
            sorted(tf.relation_to_id.keys()),
        )

    return tf
