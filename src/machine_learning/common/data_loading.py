"""
src/machine_learning/common/data_loading.py - BioCypher CSV loading shared by KGE and R-GCN.

Provides low-level helpers to read BioCypher header+part CSV pairs and produce
a (head, relation, tail) triple array consumable by both PyKEEN TriplesFactory
construction (KGE track) and PyG heterogeneous graph building (R-GCN track).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

CSV_SEP = '\t'

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_BIOCYPHER_ROOT = _PROJECT_ROOT / 'data' / 'biocypher'


# ── Low-level CSV reader ──────────────────────────────────────────────────────


def _read_type(
    bc_dir: Path,
    type_name: str,
    sep: str = CSV_SEP,
) -> pd.DataFrame | None:
    """
    Read all part files for *type_name*; returns None if no files were found.
    """
    header_path = bc_dir / f'{type_name}-header.csv'
    if not header_path.exists():
        logger.warning('Header file not found: %s', header_path)
        return None

    columns: list[str] = (
        pd.read_csv(header_path, sep=sep, header=None, dtype=str).iloc[0].tolist()
    )

    part_files = sorted(bc_dir.glob(f'{type_name}-part*.csv'))
    if not part_files:
        logger.warning(
            "No data part files found for type '%s' in %s", type_name, bc_dir
        )
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
            logger.error('Could not read %s: %s', part, exc)

    if not chunks:
        return None

    result = pd.concat(chunks, ignore_index=True)
    logger.debug(
        '  %s - %d rows from %d part file(s)', type_name, len(result), len(part_files)
    )
    return result


# ── Public API ────────────────────────────────────────────────────────────────


def resolve_biocypher_dir(
    data_dir: str | Path | None = None,
    species: str = 'human',
) -> Path:
    """
    Return an explicit BioCypher dir, or the newest data/biocypher/<species>_* match.

    Args:
        data_dir: Explicit path. When given, returned as-is.
        species: Species prefix used to find the newest auto-detected build
            (e.g. 'human' matches data/biocypher/human_20241201_120000/).

    Returns:
        Resolved Path to the BioCypher output directory.

    Raises:
        FileNotFoundError: If auto-detection finds no candidates.
    """
    if data_dir is not None:
        return Path(data_dir)

    candidates = sorted(
        _BIOCYPHER_ROOT.glob(f'{species}_*'),
        key=lambda p: (p.stat().st_mtime, p.name),
        reverse=True,
    )
    if candidates:
        return candidates[0]

    raise FileNotFoundError(
        f'No BioCypher output directory found under {_BIOCYPHER_ROOT} for {species!r}.\n'
        f'Run: uv run python run_kg.py --species {species}'
    )


def discover_edge_types(bc_dir: Path, sep: str = CSV_SEP) -> list[str]:
    """
    Return all edge type names that have both a header file and at least one part file.

    Edge headers are identified by the presence of both ':START_ID' and ':END_ID'
    columns (node headers have ':ID' instead).

    Args:
        bc_dir: BioCypher output directory.
        sep: CSV field separator.

    Returns:
        Sorted list of edge type name strings.
    """
    edge_types: list[str] = []

    for header_file in sorted(bc_dir.glob('*-header.csv')):
        type_name = header_file.stem.replace('-header', '')
        try:
            cols: list[str] = (
                pd.read_csv(header_file, sep=sep, header=None, dtype=str)
                .iloc[0]
                .tolist()
            )
        except Exception:
            continue

        if ':START_ID' in cols and ':END_ID' in cols:
            if list(bc_dir.glob(f'{type_name}-part*.csv')):
                edge_types.append(type_name)

    return edge_types


def load_triples(
    bc_dir: str | Path | None,
    include_relations: set[str] | None = None,
    sep: str = CSV_SEP,
) -> np.ndarray:
    """
    Load all edges from BioCypher CSVs as (head, relation, tail) string triples.

    include_relations=None (default) loads all relation types.

    Args:
        bc_dir: BioCypher output directory (None → auto-detect).
        include_relations: Relation type names to include. None loads all.
        sep: CSV field separator.

    Returns:
        np.ndarray of shape (N, 3): [head_id, relation_type, tail_id].

    Raises:
        FileNotFoundError: If bc_dir does not exist.
        RuntimeError: If no edge files or no triples were loaded.
    """
    bc_dir = resolve_biocypher_dir(bc_dir)
    if not bc_dir.exists():
        raise FileNotFoundError(
            f'BioCypher output directory not found: {bc_dir}\n'
            'Run the build pipeline first: uv run python run_kg.py'
        )

    edge_types = discover_edge_types(bc_dir, sep=sep)
    if not edge_types:
        raise RuntimeError(
            f'No edge files found in {bc_dir}.\n'
            'Check that BioCypher wrote *-part*.csv files alongside the headers.'
        )

    logger.info('Found %d edge type(s) in %s', len(edge_types), bc_dir)

    all_arrays: list[np.ndarray] = []

    for etype in edge_types:
        if include_relations is not None and etype not in include_relations:
            logger.debug('  Skipping relation type: %s', etype)
            continue

        df = _read_type(bc_dir, etype, sep=sep)
        if df is None or df.empty:
            logger.warning('  No data for edge type: %s - skipping', etype)
            continue

        if ':TYPE' not in df.columns:
            df[':TYPE'] = etype

        triples = df[[':START_ID', ':TYPE', ':END_ID']].dropna().values.astype(str)

        # BioCypher wraps :TYPE values in single quotes - strip them.
        triples[:, 1] = np.char.strip(triples[:, 1], "'")

        all_arrays.append(triples)
        logger.info('  %s: %d triples', etype, len(triples))

    if not all_arrays:
        raise RuntimeError(
            f'No triples were loaded. include_relations filter was: {include_relations}'
        )

    combined = np.concatenate(all_arrays, axis=0)
    logger.info('Total triples: %d', len(combined))
    return combined
