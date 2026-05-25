"""
src/knowledge_graph/layers/mouse_gene_to_phenotype_export.py - mouse gene-to-phenotype helpers.

Provides shared functions for reading MGI phenotype data and collapsing
mouse gene→MP term associations to MP top-level terms.  Called by
human_gene_to_phenotype_export.py for the mouse branch.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from src.knowledge_graph.gene_ids import validate_required_columns
from src.knowledge_graph.layers._mp_ontology import collapse_to_top

logger = logging.getLogger(__name__)


def _load_mp_ontology(data_dir: Path):
    """
    Load the MP ontology from mp.obo in data_dir using pronto.

    Args:
        data_dir: Directory containing mp.obo.

    Returns:
        A pronto.Ontology instance.

    Raises:
        ImportError: If pronto is not installed.
        FileNotFoundError: If mp.obo is absent.
    """
    try:
        from pronto import Ontology
    except ImportError as exc:
        raise ImportError(
            'pronto is required for MP ontology parsing.\nInstall with: uv add pronto'
        ) from exc

    mp_path = data_dir / 'mp.obo'
    if not mp_path.exists():
        raise FileNotFoundError(
            f'mp.obo not found: {mp_path}\n'
            'Download from http://www.informatics.jax.org/downloads/reports/mp.obo'
        )
    logger.debug('Loading MP ontology from %s', mp_path)
    return Ontology(str(mp_path))


def _build_mp_to_top_edges(
    mp_terms: pd.Series,
    ont,
    top_set: set[str],
    root_id: str,
) -> pd.DataFrame:
    """
    Map a series of MP term IDs to their top-level ancestor IDs.

    Args:
        mp_terms: Series of MP term ID strings.
        ont: pronto.Ontology instance.
        top_set: Set of top-level MP term IDs.
        root_id: Root MP term ID to exclude from results.

    Returns:
        DataFrame with columns: mp_id, mp_top_id.
    """
    rows = []
    for tid in sorted(t for t in mp_terms.dropna().unique() if isinstance(t, str)):
        for top in collapse_to_top(ont, tid, top_set, root_id):
            rows.append({'mp_id': tid, 'mp_top_id': top})
    return (
        pd.DataFrame(rows, columns=['mp_id', 'mp_top_id'])
        if rows
        else pd.DataFrame(columns=['mp_id', 'mp_top_id'])
    )


def _read_mgi_phenogenomp(path: str | Path) -> pd.DataFrame:
    """
    Read the MGI PhenoGenoMP.rpt file.

    Args:
        path: Path to MGI_PhenoGenoMP.rpt.

    Returns:
        DataFrame with columns: genotype, allele, background, mp_id, pmid,
        marker_mgi, genotype_mgi.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    cols = [
        'genotype',
        'allele',
        'background',
        'mp_id',
        'pmid',
        'marker_mgi',
        'genotype_mgi',
    ]
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f'MGI_PhenoGenoMP.rpt not found: {path}\n'
            'Download from https://www.informatics.jax.org/downloads/reports/'
        )
    return pd.read_csv(path, sep='\t', header=None, names=cols, dtype=str)


def _build_mouse_gene_mp_top_edges(
    mouse_gene_mp: pd.DataFrame,
    mp_to_top: pd.DataFrame,
) -> pd.DataFrame:
    """
    Join mouse gene→MP edges with an MP→top mapping to produce gene→MP top edges.

    Args:
        mouse_gene_mp: DataFrame with columns gene_id, gene_symbol, mp_id.
        mp_to_top: DataFrame with columns mp_id, mp_top_id.

    Returns:
        Deduplicated DataFrame with columns gene_id, gene_symbol, mp_top_id.
    """
    validate_required_columns(
        mouse_gene_mp,
        {'gene_id', 'gene_symbol', 'mp_id'},
        'mouse phenotype export input',
    )
    validate_required_columns(mp_to_top, {'mp_id', 'mp_top_id'}, 'mouse MP top mapping')
    return (
        mouse_gene_mp[['gene_id', 'gene_symbol', 'mp_id']]
        .drop_duplicates()
        .merge(mp_to_top, on='mp_id', how='inner')[
            ['gene_id', 'gene_symbol', 'mp_top_id']
        ]
        .drop_duplicates()
        .sort_values(['gene_id', 'mp_top_id'])
        .reset_index(drop=True)
    )
