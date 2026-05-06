from __future__ import annotations

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
    """Read all part files for *type_name* and return a single DataFrame.

    Parameters
    ----------
    bc_dir    : directory containing BioCypher output (e.g. biocypher_out/human/)
    type_name : BioCypher type name, e.g. 'HumanGeneHasMpTopTerm'
    sep       : CSV field separator (default: tab)

    Returns
    -------
    DataFrame with the header columns, or None if no files were found.
    """
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

    Parameters
    ----------
    bc_dir             : biocypher_out/human/ directory
    include_relations  : if given, only these relation types are loaded
                         (e.g. {"HumanGeneHasMpTopTerm"}).
                         None (default) = load everything — recommended, because
                         the richer context improves all entity embeddings.
    sep                : CSV separator (default: tab, from biocypher_config.yaml)

    Returns
    -------
    np.ndarray of shape (N, 3) with dtype str.
    Columns: [head_id, relation_type, tail_id]
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


def build_triples_factory(
    bc_dir: str | Path,
    include_relations: Optional[set[str]] = None,
    create_inverse_triples: bool = False,
    sep: str = CSV_SEP,
) -> TriplesFactory:
    """Build a PyKEEN TriplesFactory from a BioCypher output directory.

    Parameters
    ----------
    bc_dir                   : biocypher_out/human/
    include_relations        : relation whitelist (None = all)
    create_inverse_triples   : if True, add (tail, relation_inverse, head) for
                               every triple — can improve tail prediction quality
                               at the cost of doubling the triple count.
    sep                      : CSV separator

    Returns
    -------
    TriplesFactory ready to be split and fed to PyKEEN models.
    """
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

    # Sanity check: warn if the primary prediction relation is missing
    if GENE_PHENOTYPE_RELATION not in tf.relation_to_id:
        logger.warning(
            "Primary relation '%s' not found in the loaded data!\n"
            "Available relations: %s",
            GENE_PHENOTYPE_RELATION,
            sorted(tf.relation_to_id.keys()),
        )

    return tf