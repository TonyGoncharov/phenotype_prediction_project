from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from pykeen.triples import TriplesFactory

from .config import CSV_SEP, GENE_PHENOTYPE_RELATION

logger = logging.getLogger(__name__)


# ── Low-level helpers ────────────────────────────────────────────────────────


def _read_type(
    bc_dir: Path,
    type_name: str,
    sep: str = CSV_SEP,
) -> Optional[pd.DataFrame]:
    """Read all part files for *type_name*; returns None if no files were found."""
    header_path = bc_dir / f"{type_name}-header.csv"
    if not header_path.exists():
        logger.warning("Header file not found: %s", header_path)
        return None

    # Header file contains a single row with column names
    columns: list[str] = (
        pd.read_csv(header_path, sep=sep, header=None, dtype=str)
        .iloc[0]
        .tolist()
    )

    part_files = sorted(bc_dir.glob(f"{type_name}-part*.csv"))
    if not part_files:
        logger.warning("No data part files found for type '%s' in %s", type_name, bc_dir)
        return None

    chunks: list[pd.DataFrame] = []
    for part in part_files:
        try:
            df = pd.read_csv(
                part,
                sep=sep,
                header=None,
                names=columns,
                dtype=str,
                low_memory=False,
            )
            chunks.append(df)
        except Exception as exc:
            logger.error("Could not read %s: %s", part, exc)

    if not chunks:
        return None

    result = pd.concat(chunks, ignore_index=True)
    logger.debug("  %s — %d rows from %d part file(s)", type_name, len(result), len(part_files))
    return result


# ── Public API ────────────────────────────────────────────────────────────────


def discover_edge_types(bc_dir: Path, sep: str = CSV_SEP) -> list[str]:
    """Return all *edge* type names that have both a header and at least one part file.

    Edge headers are identified by the presence of both ':START_ID' and ':END_ID'
    columns (as opposed to node headers which have ':ID').
    """
    edge_types: list[str] = []

    for header_file in sorted(bc_dir.glob("*-header.csv")):
        type_name = header_file.stem.replace("-header", "")
        try:
            cols: list[str] = (
                pd.read_csv(header_file, sep=sep, header=None, dtype=str)
                .iloc[0]
                .tolist()
            )
        except Exception:
            continue

        if ":START_ID" in cols and ":END_ID" in cols:
            if list(bc_dir.glob(f"{type_name}-part*.csv")):
                edge_types.append(type_name)

    return edge_types


def load_triples(
    bc_dir: str | Path,
    include_relations: Optional[set[str]] = None,
    sep: str = CSV_SEP,
) -> np.ndarray:
    """Load all edges from BioCypher CSVs as (head, relation, tail) string triples.

    include_relations=None (default) loads all relation types — recommended
    because the richer context improves all entity embeddings.

    Returns np.ndarray of shape (N, 3): [head_id, relation_type, tail_id].
    """
    bc_dir = Path(bc_dir)
    if not bc_dir.exists():
        raise FileNotFoundError(
            f"BioCypher output directory not found: {bc_dir}\n"
            "Run the build pipeline first:  python run.py --species human"
        )

    edge_types = discover_edge_types(bc_dir, sep=sep)
    if not edge_types:
        raise RuntimeError(
            f"No edge files found in {bc_dir}.\n"
            "Check that BioCypher wrote *-part*.csv files alongside the headers."
        )

    logger.info("Found %d edge type(s) in %s", len(edge_types), bc_dir)

    all_arrays: list[np.ndarray] = []

    for etype in edge_types:
        if include_relations is not None and etype not in include_relations:
            logger.debug("  Skipping relation type: %s", etype)
            continue

        df = _read_type(bc_dir, etype, sep=sep)
        if df is None or df.empty:
            logger.warning("  No data for edge type: %s — skipping", etype)
            continue

        # If :TYPE column is absent (shouldn't happen), fall back to the file name
        if ":TYPE" not in df.columns:
            df[":TYPE"] = etype

        triples = (
            df[[":START_ID", ":TYPE", ":END_ID"]]
            .dropna()
            .values
            .astype(str)
        )

        # BioCypher wraps :TYPE values in single quotes, e.g. 'HumanGeneHasMpTopTerm'
        # Strip them so relation names match what we expect in config.py
        triples[:, 1] = np.char.strip(triples[:, 1], "'")

        all_arrays.append(triples)
        logger.info("  %-42s  %7d triples", etype, len(triples))

    if not all_arrays:
        raise RuntimeError(
            "No triples were loaded.  "
            f"include_relations filter was: {include_relations}"
        )

    combined = np.concatenate(all_arrays, axis=0)
    logger.info("Total triples: %d", len(combined))
    return combined


def _tf_cache_key(
    include_relations: Optional[set[str]],
    create_inverse_triples: bool,
    source_files: frozenset[str] = frozenset(),
) -> str:
    """Short deterministic key encoding the loader parameters and source file set.

    *source_files* is the set of CSV filenames in the source directory — including
    it ensures the cache is invalidated when new files are added, not only when
    existing files are modified.
    """
    payload = json.dumps(
        {
            "rels":  sorted(include_relations) if include_relations else None,
            "inv":   create_inverse_triples,
            "files": sorted(source_files),
        },
        sort_keys=True,
    )
    return hashlib.md5(payload.encode()).hexdigest()[:10]


def _cache_is_fresh(cache_path: Path, source_dir: Path) -> bool:
    """True if *cache_path* exists and is newer than every CSV in *source_dir*."""
    if not cache_path.exists():
        return False
    cache_mtime = cache_path.stat().st_mtime
    return all(f.stat().st_mtime <= cache_mtime for f in source_dir.glob("*.csv"))


def build_triples_factory(
    bc_dir: str | Path,
    include_relations: Optional[set[str]] = None,
    create_inverse_triples: bool = False,
    sep: str = CSV_SEP,
    cache_dir: Optional[str | Path] = None,
) -> TriplesFactory:
    """Build a PyKEEN TriplesFactory from a BioCypher output directory.

    create_inverse_triples adds (tail, relation_inverse, head) for every triple —
    can improve tail prediction quality at the cost of doubling the triple count.

    cache_dir persists the TriplesFactory as a binary and reloads it on subsequent
    calls; invalidated automatically when any source CSV is newer than the cache or
    when the set of CSV files in bc_dir changes (files added or removed).
    """
    bc_dir = Path(bc_dir)

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        source_files = frozenset(f.name for f in bc_dir.glob("*.csv"))
        key = _tf_cache_key(include_relations, create_inverse_triples, source_files)
        cache_path = cache_dir / f"tf_{key}"
        if _cache_is_fresh(cache_path, bc_dir):
            logger.info("Loading TriplesFactory from cache: %s", cache_path)
            tf = TriplesFactory.from_path_binary(cache_path)
            logger.info(
                "TriplesFactory ready (cached) — %d entities | %d relations | %d triples",
                tf.num_entities, tf.num_relations, tf.num_triples,
            )
            return tf
        logger.info("Cache miss (key=%s) — building from CSV files", key)

    triples = load_triples(bc_dir, include_relations=include_relations, sep=sep)

    tf = TriplesFactory.from_labeled_triples(
        triples=triples,
        create_inverse_triples=create_inverse_triples,
    )

    logger.info(
        "TriplesFactory ready — %d entities | %d relations | %d triples",
        tf.num_entities,
        tf.num_relations,
        tf.num_triples,
    )

    if cache_dir is not None:
        tf.to_path_binary(cache_path)
        logger.info("TriplesFactory cached to: %s", cache_path)

    # Sanity check: warn if the primary prediction relation is missing
    if GENE_PHENOTYPE_RELATION not in tf.relation_to_id:
        logger.warning(
            "Primary relation '%s' not found in the loaded data!\n"
            "Available relations: %s",
            GENE_PHENOTYPE_RELATION,
            sorted(tf.relation_to_id.keys()),
        )

    return tf