"""
src/knowledge_graph/layers/human_gene_expression_export.py - GTEx expression export pipeline.

Reads GTEx median-TPM GCT files, applies a TPM threshold, and writes
BioCypher-ready node/edge TSVs for the human gene-expression layer.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd
from src.knowledge_graph.gene_ids import human_gene_id

logger = logging.getLogger(__name__)

DEFAULT_TPM_THRESHOLD = 5.0

_GCT_ID_COL = 'Name'
_GCT_SYM_COL = 'Description'


def _gtex_tissue_id(tissue_name: str) -> str:
    """
    Convert a GTEx tissue name to a stable GTEX:<slug> identifier.

    Args:
        tissue_name: Raw GTEx tissue name string.

    Returns:
        Normalised identifier of the form 'GTEX:<slug>'.
    """
    return 'GTEX:' + re.sub(r'[^a-zA-Z0-9]+', '_', tissue_name).strip('_')


def read_gtex_gct(path: str | Path) -> pd.DataFrame:
    """
    Parse a GTEx GCT median-TPM file into long format.

    Args:
        path: Path to the GCT or GCT.gz file.

    Returns:
        Long-format DataFrame with columns: gene_symbol, tissue_name, median_tpm.

    Raises:
        FileNotFoundError: If the GCT file does not exist.
        ValueError: If expected GCT columns are absent.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f'GTEx GCT file not found: {path}\n'
            'Download from https://gtexportal.org/home/datasets\n'
            '  (GTEx_Analysis_v*_gene_median_tpm.gct.gz)'
        )

    logger.debug('Parsing GCT file: %s', path)

    df = pd.read_csv(path, sep='\t', skiprows=2, dtype=str)

    if _GCT_ID_COL not in df.columns or _GCT_SYM_COL not in df.columns:
        raise ValueError(
            f"Unexpected GCT format - expected columns '{_GCT_ID_COL}' and "
            f"'{_GCT_SYM_COL}', got: {list(df.columns[:4])} ..."
        )

    tissue_cols = [c for c in df.columns if c not in (_GCT_ID_COL, _GCT_SYM_COL)]
    if not tissue_cols:
        raise ValueError('No tissue columns found in GCT file.')
    logger.debug('GCT: %d gene rows, %d tissue columns', len(df), len(tissue_cols))

    df_long = df.melt(
        id_vars=[_GCT_SYM_COL],
        value_vars=tissue_cols,
        var_name='tissue_name',
        value_name='median_tpm',
    ).rename(columns={_GCT_SYM_COL: 'gene_symbol'})

    before = len(df_long)
    df_long['median_tpm'] = pd.to_numeric(df_long['median_tpm'], errors='coerce')
    df_long = df_long.dropna(subset=['gene_symbol', 'tissue_name', 'median_tpm'])
    if dropped := before - len(df_long):
        logger.warning('Dropped %d rows with missing/non-numeric TPM values', dropped)

    return df_long[['gene_symbol', 'tissue_name', 'median_tpm']]


def run_expression_pipeline(
    gtex_path: str | Path,
    out_dir: str | Path = './out',
    tpm_threshold: float = DEFAULT_TPM_THRESHOLD,
) -> dict[str, pd.DataFrame]:
    """
    Run the human GTEx expression export pipeline.

    Args:
        gtex_path: Path to the GTEx GCT file.
        out_dir: Output directory for TSV files.
        tpm_threshold: Minimum median TPM for a gene-tissue edge to be included.

    Returns:
        Dict with keys 'edges' and 'tissues' mapping to the output DataFrames.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Starting GTEx expression pipeline')
    logger.info('GTEx file: %s', gtex_path)
    logger.info('Output dir: %s', out_dir)
    logger.info('TPM threshold: %.1f', tpm_threshold)

    expr = read_gtex_gct(gtex_path)
    all_symbols = set(expr['gene_symbol'].dropna().unique())
    logger.info('Genes in GTEx: %d', len(all_symbols))
    logger.info('Tissues in GTEx: %d', expr['tissue_name'].nunique())

    expressed = expr[expr['median_tpm'] >= tpm_threshold].copy()
    expressed['tissue_id'] = expressed['tissue_name'].map(_gtex_tissue_id)
    # MVP convention: HGNC:<symbol> is not a canonical HGNC identifier.
    # TODO: migrate human gene nodes to true HGNC IDs when available.
    expressed['gene_id'] = expressed['gene_symbol'].map(human_gene_id)

    edge_df = expressed.groupby(
        ['gene_id', 'gene_symbol', 'tissue_id', 'tissue_name'], as_index=False
    )['median_tpm'].max()
    edge_df.to_csv(
        out_dir / 'edge_human_gene_expressed_in_tissue.tsv', sep='\t', index=False
    )
    logger.info(
        'Written: edge_human_gene_expressed_in_tissue.tsv (%d rows)', len(edge_df)
    )
    if d := edge_df.duplicated(subset=['gene_id', 'tissue_id']).sum():
        logger.warning('Duplicate gene→tissue rows in output: %d', d)

    tissue_df = (
        expressed[['tissue_id', 'tissue_name']]
        .drop_duplicates()
        .sort_values('tissue_name')
        .reset_index(drop=True)
    )
    tissue_df.to_csv(out_dir / 'node_tissue_gtex.tsv', sep='\t', index=False)
    logger.info('Written: node_tissue_gtex.tsv (%d rows)', len(tissue_df))

    unexpressed = sorted(all_symbols - set(edge_df['gene_symbol'].dropna()))
    pd.DataFrame({'gene_symbol': unexpressed}).to_csv(
        out_dir / 'qc_expression_unexpressed_genes.tsv', sep='\t', index=False
    )
    logger.info('Gene→tissue edges: %d', len(edge_df))
    logger.info('Unique genes expressed: %d', edge_df['gene_id'].nunique())
    logger.info('Unique tissues: %d', edge_df['tissue_id'].nunique())
    logger.info('Genes below threshold: %d', len(unexpressed))

    return {'edges': edge_df, 'tissues': tissue_df}
