"""
src/knowledge_graph/layers/ppi_export.py - BioGRID protein-protein interaction export pipeline.

Reads a BioGRID tab3-format file, filters to physical interactions for a single
organism, and writes protein node, gene-encodes-protein edge, and PPI edge TSVs
ready for BioCypher import.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from src.knowledge_graph.gene_ids import (
    human_gene_id,
    mouse_gene_id,
    validate_required_columns,
)

logger = logging.getLogger(__name__)

PROTEIN_ID_PREFIX = 'PROTEIN:'


@dataclass(frozen=True)
class OrganismConfig:
    """
    Configuration bundle for a single organism's PPI export.
    """

    tax_id: str
    organism: str
    file_prefix: str
    encodes_label: str
    ppi_label: str


HUMAN = OrganismConfig(
    tax_id='9606',
    organism='Homo sapiens',
    file_prefix='human',
    encodes_label='human gene encodes protein',
    ppi_label='human protein interacts with',
)

MOUSE = OrganismConfig(
    tax_id='10090',
    organism='Mus musculus',
    file_prefix='mouse',
    encodes_label='mouse gene encodes protein',
    ppi_label='mouse protein interacts with',
)


_BIOGRID_FIELD_ALIASES: dict[str, tuple[str, str]] = {
    'symbol_a': ('Official Symbol Interactor A', 'Gene Name Interactor A'),
    'symbol_b': ('Official Symbol Interactor B', 'Gene Name Interactor B'),
    'experimental_system': ('Experimental System', 'Experimental System'),
    'experimental_system_type': (
        'Experimental System Type',
        'Experimental System Type',
    ),
    'pubmed_id': ('Publication Source', 'NCBI Publication ID'),
    'tax_id_a': ('Organism ID Interactor A', 'Organism Interactor A'),
    'tax_id_b': ('Organism ID Interactor B', 'Organism Interactor B'),
    'throughput': ('Throughput', 'Throughput'),
    'source_db': ('Source Database', 'Source Database'),
}


def _detect_column_map(columns: list[str]) -> dict[str, str]:
    """
    Map BioGRID raw column names to internal field names.

    Supports both v4 and v5 tab3 column naming conventions.

    Args:
        columns: List of column names from the raw BioGRID file.

    Returns:
        Dict mapping raw column name → internal field name.

    Raises:
        ValueError: If required fields cannot be mapped to any known column name.
    """
    col_set = set(columns)

    for gen_idx in (0, 1):
        rename: dict[str, str] = {}
        for internal, aliases in _BIOGRID_FIELD_ALIASES.items():
            raw = aliases[gen_idx]
            if raw in col_set:
                rename[raw] = internal
        if len(rename) == len(_BIOGRID_FIELD_ALIASES):
            logger.debug(
                'BioGRID column format detected: %s', 'v4' if gen_idx == 0 else 'v5'
            )
            return rename

    missing_fields = {
        f: aliases
        for f, aliases in _BIOGRID_FIELD_ALIASES.items()
        if aliases[0] not in col_set and aliases[1] not in col_set
    }
    raise ValueError(
        f'Cannot map {len(missing_fields)} required BioGRID fields to any known column name.\n\n'
        'Missing fields and their expected names:\n'
        + '\n'.join(
            f"  {field}:\n    v4 → '{v4}'\n    v5 → '{v5}'"
            for field, (v4, v5) in _BIOGRID_FIELD_ALIASES.items()
            if field in missing_fields
        )
        + f'\n\nActual columns in file ({len(columns)}):\n  '
        + '\n  '.join(columns[:40])
        + ('\n  …' if len(columns) > 40 else '')
        + '\n\nMake sure you are using a BioGRID tab3-format file.'
    )


def read_biogrid(path: str | Path) -> pd.DataFrame:
    """
    Read a BioGRID tab3 file and normalise its column names.

    Args:
        path: Path to the BioGRID tab3 file.

    Returns:
        DataFrame with internal field names as columns.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the column mapping cannot be determined.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f'BioGRID file not found: {path}\n'
            'Download from https://downloads.thebiogrid.org/BioGRID/'
        )
    logger.info('Reading BioGRID from %s', path)

    df = pd.read_csv(path, sep='\t', dtype=str, low_memory=False)
    df.columns = [c.lstrip('#') for c in df.columns]
    logger.debug('BioGRID raw rows: %d, columns: %d', len(df), len(df.columns))

    rename = _detect_column_map(list(df.columns))
    df = df[list(rename)].rename(columns=rename)

    for col in ('tax_id_a', 'tax_id_b'):
        df[col] = (
            df[col].str.extract(r'\((\d+)\)', expand=False).fillna(df[col]).str.strip()
        )

    logger.debug('BioGRID rows after column normalisation: %d', len(df))
    return df


def filter_physical(df: pd.DataFrame, cfg: OrganismConfig) -> pd.DataFrame:
    """
    Filter BioGRID rows to same-organism physical interactions.

    Args:
        df: Normalised BioGRID DataFrame.
        cfg: OrganismConfig for the target species.

    Returns:
        Filtered DataFrame containing only physical intra-species interactions.
    """
    before = len(df)
    df = df[(df['tax_id_a'] == cfg.tax_id) & (df['tax_id_b'] == cfg.tax_id)].copy()
    logger.debug('After %s-only filter: %d / %d rows', cfg.organism, len(df), before)

    df = df[df['experimental_system_type'].str.lower() == 'physical'].copy()
    logger.debug('After physical-only filter: %d rows', len(df))

    df = df.dropna(subset=['symbol_a', 'symbol_b'])
    df = df[df['symbol_a'] != df['symbol_b']]
    logger.debug('After drop-NA / self-interactions: %d rows', len(df))
    return df


def build_protein_nodes(ppi: pd.DataFrame, cfg: OrganismConfig) -> pd.DataFrame:
    """
    Build a protein node DataFrame from unique interactor symbols.

    Args:
        ppi: Filtered BioGRID DataFrame.
        cfg: OrganismConfig for the target species.

    Returns:
        DataFrame with columns: protein_id, gene_symbol, organism.
    """
    validate_required_columns(ppi, {'symbol_a', 'symbol_b'}, 'BioGRID PPI input')
    symbols = pd.concat([ppi['symbol_a'], ppi['symbol_b']]).dropna().unique()
    logger.debug('Unique protein symbols (%s): %d', cfg.organism, len(symbols))
    return pd.DataFrame(
        {
            'protein_id': [f'{PROTEIN_ID_PREFIX}{s}' for s in symbols],
            'gene_symbol': list(symbols),
            'organism': cfg.organism,
        }
    )


def build_gene_encodes_protein(
    protein_nodes: pd.DataFrame,
    cfg: OrganismConfig,
    gene_id_map: dict[str, str] | None = None,
) -> pd.DataFrame:
    """
    Build gene-encodes-protein edges from protein node records.

    Args:
        protein_nodes: DataFrame with gene_symbol and protein_id columns.
        cfg: OrganismConfig for the target species.
        gene_id_map: Optional pre-built symbol → gene_id mapping.

    Returns:
        DataFrame with columns: gene_id, gene_symbol, protein_id.
    """
    validate_required_columns(
        protein_nodes, {'gene_symbol', 'protein_id'}, 'protein node input'
    )
    df = protein_nodes[['gene_symbol', 'protein_id']].copy()

    if gene_id_map is not None:
        df['gene_id'] = df['gene_symbol'].map(gene_id_map)
        n_missing = df['gene_id'].isna().sum()
        if n_missing:
            logger.warning(
                '%d symbol(s) not in gene_id_map - omitted from encodes edges '
                '(see QC file).',
                n_missing,
            )
        df = df.dropna(subset=['gene_id'])
    else:
        # Human MVP IDs use HGNC:<symbol>, which is not a canonical HGNC ID.
        # TODO: migrate human to true HGNC IDs and mouse to true MGI IDs later.
        formatter = human_gene_id if cfg.file_prefix == 'human' else mouse_gene_id
        df['gene_id'] = df['gene_symbol'].map(formatter)

    return df[['gene_id', 'gene_symbol', 'protein_id']].reset_index(drop=True)


def build_ppi_edges(
    ppi: pd.DataFrame,
    cfg: OrganismConfig,
    min_publications: int = 1,
) -> pd.DataFrame:
    """
    Aggregate BioGRID rows into deduplicated PPI edges with metadata.

    Args:
        ppi: Filtered BioGRID DataFrame.
        cfg: OrganismConfig for the target species.
        min_publications: Minimum distinct PubMed IDs required to retain an edge.

    Returns:
        DataFrame with columns: protein_id_a, protein_id_b, experimental_systems,
        experimental_system_types, throughput, pubmed_ids, source_db.
    """
    validate_required_columns(
        ppi,
        {
            'symbol_a',
            'symbol_b',
            'experimental_system',
            'experimental_system_type',
            'throughput',
            'pubmed_id',
            'source_db',
        },
        'BioGRID PPI input',
    )
    df = ppi.copy()
    df['protein_id_a'] = PROTEIN_ID_PREFIX + df['symbol_a']
    df['protein_id_b'] = PROTEIN_ID_PREFIX + df['symbol_b']

    mask = df['protein_id_a'] > df['protein_id_b']
    df.loc[mask, ['protein_id_a', 'protein_id_b']] = df.loc[
        mask, ['protein_id_b', 'protein_id_a']
    ].values

    def _sorted_unique(series: pd.Series) -> list[str]:
        """
        Return sorted unique non-null values from a Series.
        """
        return sorted(series.dropna().unique().tolist())

    def _throughput(series: pd.Series) -> str:
        """
        Classify throughput as Low if any row is low-throughput, else High.
        """
        vals = series.dropna().str.lower()
        return (
            'Low Throughput' if (vals == 'low throughput').any() else 'High Throughput'
        )

    agg = (
        df.groupby(['protein_id_a', 'protein_id_b'])
        .agg(
            experimental_systems=('experimental_system', _sorted_unique),
            experimental_system_types=('experimental_system_type', _sorted_unique),
            throughput=('throughput', _throughput),
            pubmed_ids=('pubmed_id', _sorted_unique),
            source_db=('source_db', 'first'),
        )
        .reset_index()
    )
    logger.debug('Unique PPI edges after aggregation (%s): %d', cfg.organism, len(agg))

    if min_publications >= 1:
        before = len(agg)
        agg = agg[agg['pubmed_ids'].apply(len) >= min_publications].reset_index(
            drop=True
        )
        logger.info(
            'min_publications=%d filter: %d → %d edges (removed %d, %.1f%%)',
            min_publications,
            before,
            len(agg),
            before - len(agg),
            (before - len(agg)) / before * 100,
        )

    return agg


def run_ppi_pipeline(
    biogrid_path: str | Path,
    cfg: OrganismConfig,
    out_dir: str | Path = './out',
    min_publications: int = 1,
    gene_id_map: dict[str, str] | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Run the BioGRID PPI export pipeline for a single organism.

    Args:
        biogrid_path: Path to the BioGRID tab3 file.
        cfg: OrganismConfig specifying the target organism.
        out_dir: Output directory for TSV files.
        min_publications: Minimum distinct PubMed IDs per PPI pair.
        gene_id_map: Optional symbol → gene_id map; derived from cfg if None.

    Returns:
        Dict with keys: protein_nodes, gene_encodes_protein, ppi_edges.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    p = cfg.file_prefix

    logger.info('=' * 70)
    logger.info('PPI export pipeline (%s)', cfg.organism)
    logger.info('BioGRID path: %s', biogrid_path)
    logger.info('Output dir: %s', out_dir)
    logger.info('min_publications: %d', min_publications)
    logger.info('=' * 70)

    raw = read_biogrid(biogrid_path)
    filtered = filter_physical(raw, cfg)

    if filtered.empty:
        logger.warning(
            'No %s physical interactions found - output files will be empty.',
            cfg.organism,
        )

    protein_nodes = build_protein_nodes(filtered, cfg)
    gene_encodes_protein = build_gene_encodes_protein(protein_nodes, cfg, gene_id_map)
    ppi_edges = build_ppi_edges(filtered, cfg, min_publications)

    protein_nodes.to_csv(out_dir / f'node_{p}_protein.tsv', sep='\t', index=False)
    logger.info('Written: node_%s_protein.tsv (%d rows)', p, len(protein_nodes))

    gene_encodes_protein.to_csv(
        out_dir / f'edge_{p}_gene_encodes_protein.tsv', sep='\t', index=False
    )
    logger.info(
        'Written: edge_%s_gene_encodes_protein.tsv (%d rows)',
        p,
        len(gene_encodes_protein),
    )
    _dupes = gene_encodes_protein.duplicated(subset=['gene_id', 'protein_id']).sum()
    if _dupes:
        logger.warning('Duplicate gene→protein rows in output: %d', _dupes)

    ppi_out = ppi_edges.copy()
    for col in ('experimental_systems', 'experimental_system_types', 'pubmed_ids'):
        ppi_out[col] = ppi_out[col].apply(
            lambda v: '|'.join(v) if isinstance(v, list) else v
        )
    ppi_out.to_csv(out_dir / f'edge_{p}_ppi.tsv', sep='\t', index=False)
    logger.info('Written: edge_%s_ppi.tsv (%d rows)', p, len(ppi_out))
    _dupes = ppi_out.duplicated(subset=['protein_id_a', 'protein_id_b']).sum()
    if _dupes:
        logger.warning('Duplicate PPI rows in output: %d', _dupes)

    no_gene_id = sorted(
        set(protein_nodes['gene_symbol']) - set(gene_encodes_protein['gene_symbol'])
    )
    if no_gene_id:
        qc_path = out_dir / f'qc_{p}_ppi_symbols_no_gene_id.tsv'
        pd.DataFrame({'gene_symbol': no_gene_id}).to_csv(qc_path, sep='\t', index=False)
        logger.info(
            'QC: %d symbol(s) with no gene_id → %s', len(no_gene_id), qc_path.name
        )

    logger.info('Protein nodes: %d', len(protein_nodes))
    logger.info('Gene→Protein edges: %d', len(gene_encodes_protein))
    logger.info('Protein↔Protein edges: %d', len(ppi_edges))
    logger.info(
        'Experimental systems: %s',
        sorted({s for row in ppi_edges['experimental_systems'] for s in row}),
    )

    return {
        'protein_nodes': protein_nodes,
        'gene_encodes_protein': gene_encodes_protein,
        'ppi_edges': ppi_edges,
    }
