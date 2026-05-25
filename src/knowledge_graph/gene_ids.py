"""
src/knowledge_graph/gene_ids.py - Gene ID formatting conventions.

Defines the MVP gene ID scheme (HGNC:<symbol> for human, plain <symbol> for
mouse) and shared column-validation helpers used across all export layers.

See docs/gene_identifier_contract.md for the migration path to canonical IDs.
"""

from __future__ import annotations

import pandas as pd

HUMAN_GENE_ID_PREFIX = 'HGNC:'
HUMAN_ORGANISM = 'Homo sapiens'
MOUSE_ORGANISM = 'Mus musculus'
MOUSE_TAX_ID = '10090'


def human_gene_id(gene_symbol: str) -> str:
    """
    Return the MVP human gene ID for *gene_symbol*.

    MVP convention only: HGNC:<symbol> is not a canonical HGNC identifier.
    TODO: migrate to true HGNC IDs, e.g. HGNC:1100, when source mappings are ready.
    """
    return f'{HUMAN_GENE_ID_PREFIX}{gene_symbol}'


def mouse_gene_id(gene_symbol: str) -> str:
    """
    Return the MVP mouse gene ID for *gene_symbol*.

    MVP convention only. TODO: migrate mouse nodes to true MGI IDs later.
    """
    return gene_symbol


def add_gene_id_column(
    df: pd.DataFrame,
    species: str,
    symbol_col: str = 'gene_symbol',
    gene_id_col: str = 'gene_id',
) -> pd.DataFrame:
    """
    Add a gene_id column to *df* using the species-appropriate ID formatter.

    Args:
        df: Input DataFrame containing a gene symbol column.
        species: 'human' or 'mouse'; selects the ID formatter.
        symbol_col: Name of the column containing gene symbols.
        gene_id_col: Name of the new gene ID column to add.

    Returns:
        Copy of *df* with the new gene ID column appended.
    """
    validate_required_columns(df, {symbol_col}, 'gene ID formatting')
    out = df.copy()
    formatter = human_gene_id if species == 'human' else mouse_gene_id
    out[gene_id_col] = out[symbol_col].map(formatter)
    return out


def validate_required_columns(
    df: pd.DataFrame,
    required: set[str],
    context: str,
) -> None:
    """
    Raise ValueError if any of *required* columns are absent from *df*.

    Args:
        df: DataFrame to validate.
        required: Set of required column names.
        context: Human-readable context string for the error message.

    Raises:
        ValueError: If any required columns are missing.
    """
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(
            f'{context} missing required column(s): {missing}. Got: {list(df.columns)}'
        )
