"""
src/knowledge_graph/layers/mouse_gene_expression_export.py - GXD expression export pipeline.

Reads the MGI GXD RNA-seq archive (gxdrnaseq.tar.gz), filters for wild-type
detected expression, aggregates per gene×tissue, maps to Uberon IDs, and writes
BioCypher-ready edge TSVs for the mouse gene-expression layer.
"""

from __future__ import annotations

import gzip
import io
import logging
import tarfile
from pathlib import Path

import pandas as pd
from src.knowledge_graph.gene_ids import validate_required_columns

logger = logging.getLogger(__name__)

_USECOLS = [
    'MGI Gene ID',
    'Ensembl ID',
    'Gene Symbol',
    'Anatomical Structure',
    'Theiler Stage',
    'Mutant Allele Pair(s)',
    'Detected',
    'avg_TPM',
]


_TPM_BINS = [0, 0.5, 10, 100, float('inf')]
_TPM_LABELS = ['Absent', 'Low', 'Medium', 'High']


def _read_rpt_gz(name: str, data: bytes) -> pd.DataFrame | None:
    """
    Parse a single gzip-compressed pipe-delimited rpt file from the GXD archive.

    Args:
        name: Member filename (used for warning messages).
        data: Raw bytes of the .rpt.gz file.

    Returns:
        DataFrame with columns from _USECOLS, or None on parse failure.
    """
    try:
        with gzip.open(io.BytesIO(data), 'rb') as gz:
            raw = gz.read()
        df = pd.read_csv(
            io.BytesIO(raw),
            sep='|',
            dtype=str,
            low_memory=False,
        )
        df.columns = df.columns.str.strip()

        missing = [c for c in _USECOLS if c not in df.columns]
        if missing:
            logger.warning('Skipping %s - missing columns: %s', name, missing)
            return None

        df = df[_USECOLS].copy()
        df['avg_TPM'] = pd.to_numeric(df['avg_TPM'], errors='coerce')
        return df

    except Exception as exc:
        logger.warning('Error reading %s: %s', name, exc)
        return None


def load_tar(tar_path: str | Path) -> pd.DataFrame:
    """
    Load and concatenate all .rpt.gz experiment files from the GXD tar archive.

    Args:
        tar_path: Path to gxdrnaseq.tar.gz.

    Returns:
        Combined DataFrame of all successfully parsed experiment files.

    Raises:
        FileNotFoundError: If the archive does not exist.
        RuntimeError: If no files could be parsed.
    """
    tar_path = Path(tar_path)
    if not tar_path.exists():
        raise FileNotFoundError(
            f'GXD archive not found: {tar_path}\n'
            'Place gxdrnaseq.tar.gz in data/raw/ or pass --gxd-tar to run_kg.py.'
        )

    logger.info('Reading GXD archive: %s', tar_path)
    frames: list[pd.DataFrame] = []

    with tarfile.open(tar_path, 'r:gz') as tar:
        members = [m for m in tar.getmembers() if m.name.endswith('.rpt.gz')]
        logger.info('Experiment files in archive: %d', len(members))
        n_failed = 0
        for i, member in enumerate(members, 1):
            f = tar.extractfile(member)
            if f is None:
                n_failed += 1
                continue
            df = _read_rpt_gz(member.name, f.read())
            if df is not None:
                frames.append(df)
            else:
                n_failed += 1
            if i % 20 == 0:
                logger.debug('  Parsed %d / %d files', i, len(members))

    if n_failed:
        logger.warning(
            '%d / %d experiment files failed to parse', n_failed, len(members)
        )

    if not frames:
        raise RuntimeError('No experiment files could be parsed from the archive.')

    combined = pd.concat(frames, ignore_index=True)
    logger.info('Total rows after combining all experiments: %d', len(combined))
    return combined


def apply_filters(
    df: pd.DataFrame,
    theiler_stages: set[int] | None = None,
) -> pd.DataFrame:
    """
    Filter GXD expression rows to wild-type, detected, optionally stage-restricted.

    Args:
        df: Raw combined GXD DataFrame.
        theiler_stages: Optional set of Theiler stage numbers to retain.

    Returns:
        Filtered DataFrame.
    """
    before = len(df)

    df = df[df['Detected'].str.strip().str.upper() == 'YES']
    logger.info('After Detected=Yes filter: %d (removed %d)', len(df), before - len(df))
    before = len(df)

    wt_mask = df['Mutant Allele Pair(s)'].isna() | (
        df['Mutant Allele Pair(s)'].str.strip() == ''
    )
    df = df[wt_mask]
    logger.info('After wild-type filter: %d (removed %d)', len(df), before - len(df))
    before = len(df)

    if theiler_stages is not None:
        ts_num = pd.to_numeric(
            df['Theiler Stage'].str.extract(r'(\d+)', expand=False), errors='coerce'
        )
        df = df[ts_num.isin(theiler_stages)]
        logger.info(
            'After Theiler stage filter %s: %d (removed %d)',
            theiler_stages,
            len(df),
            before - len(df),
        )
        before = len(df)

    df = df.dropna(subset=['avg_TPM'])
    logger.info('After drop-NA TPM: %d (removed %d)', len(df), before - len(df))
    return df


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate per-experiment GXD rows to gene×tissue level statistics.

    Args:
        df: Filtered GXD DataFrame.

    Returns:
        DataFrame with one row per gene×tissue with median TPM and TPM level.
    """
    group_cols = ['MGI Gene ID', 'Ensembl ID', 'Gene Symbol', 'Anatomical Structure']
    validate_required_columns(df, set(group_cols + ['avg_TPM']), 'GXD expression input')

    agg = (
        df.groupby(group_cols, sort=False)['avg_TPM']
        .agg(n_experiments='count', median_avg_tpm='median')
        .reset_index()
    )

    agg['tpm_level'] = pd.cut(
        agg['median_avg_tpm'],
        bins=_TPM_BINS,
        labels=_TPM_LABELS,
        right=True,
        include_lowest=True,
    ).astype(str)

    agg['median_avg_tpm'] = agg['median_avg_tpm'].round(4)

    agg = agg.rename(
        columns={
            'MGI Gene ID': 'mgi_gene_id',
            'Ensembl ID': 'ensembl_id',
            'Gene Symbol': 'gene_symbol',
            'Anatomical Structure': 'anatomical_structure',
        }
    )

    logger.info('Gene × tissue pairs after aggregation: %d', len(agg))
    logger.info('Unique genes: %d', agg['mgi_gene_id'].nunique())
    logger.info('Unique tissues: %d', agg['anatomical_structure'].nunique())
    return agg


def _load_uberon_map(
    data_dir: Path,
    mapping_path: str | Path | None = None,
) -> dict[str, str]:
    """
    Load the GXD tissue name → Uberon ID mapping from TSV.

    Args:
        data_dir: Directory to search for gxd_uberon_mapping.tsv.
        mapping_path: Explicit path override; falls back to data_dir if None.

    Returns:
        Dict mapping GXD tissue name → UBERON ID for rows with a valid UBERON ID.
    """
    path = Path(mapping_path) if mapping_path else data_dir / 'gxd_uberon_mapping.tsv'
    if not path.exists():
        logger.warning(
            'gxd_uberon_mapping.tsv not found in %s - all tissues will be skipped. '
            'Run scripts/mouse/make_uberon_mapping.py first.',
            data_dir,
        )
        return {}

    df = pd.read_csv(path, sep='\t', dtype=str).fillna('')
    if not {'gxd_name', 'uberon_id'}.issubset(df.columns):
        raise ValueError(
            f"{path}: expected columns 'gxd_name' and 'uberon_id', got {list(df.columns)}"
        )

    mapped = df[df['uberon_id'].str.startswith('UBERON:', na=False)]
    logger.info('UBERON mappings loaded: %d / %d tissues mapped', len(mapped), len(df))
    return dict(zip(mapped['gxd_name'], mapped['uberon_id'], strict=True))


def run_mouse_expression_pipeline(
    tar_path: str | Path,
    data_dir: str | Path,
    out_dir: str | Path = './out',
    theiler_stages: set[int] | None = None,
    uberon_mapping_path: str | Path | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Run the mouse GXD expression export pipeline end-to-end.

    Args:
        tar_path: Path to gxdrnaseq.tar.gz.
        data_dir: Directory containing gxd_uberon_mapping.tsv.
        out_dir: Output directory for TSV files.
        theiler_stages: Optional Theiler stage filter.
        uberon_mapping_path: Explicit path to GXD→Uberon mapping TSV.

    Returns:
        Dict with keys 'edges' and 'tissues' mapping to output DataFrames.
    """
    tar_path = Path(tar_path)
    data_dir = Path(data_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Starting GXD expression pipeline (mouse)')
    logger.info('Archive: %s', tar_path)
    logger.info('Data dir: %s', data_dir)
    logger.info('Output dir: %s', out_dir)
    if theiler_stages:
        logger.info('Theiler stages filter: %s', sorted(theiler_stages))

    raw = load_tar(tar_path)
    filtered = apply_filters(raw, theiler_stages)
    summary = aggregate(filtered)

    uberon_map = _load_uberon_map(data_dir, uberon_mapping_path)
    summary['tissue_id'] = summary['anatomical_structure'].map(uberon_map)

    unmapped = sorted(
        summary[summary['tissue_id'].isna()]['anatomical_structure'].unique()
    )
    if unmapped:
        logger.warning('Skipping %d unmapped tissues (see QC file)', len(unmapped))
        pd.DataFrame({'gxd_name': unmapped}).to_csv(
            out_dir / 'qc_mouse_expression_unmapped_tissues.tsv', sep='\t', index=False
        )
    summary = summary[summary['tissue_id'].notna()].copy()

    edge_df = summary.rename(
        columns={
            'anatomical_structure': 'tissue_name',
            'median_avg_tpm': 'median_tpm',
        }
    )
    edge_df['gene_id'] = edge_df['mgi_gene_id']
    edge_df = edge_df[
        [
            'gene_id',
            'gene_symbol',
            'mgi_gene_id',
            'tissue_id',
            'tissue_name',
            'median_tpm',
            'tpm_level',
            'n_experiments',
        ]
    ]

    edge_df.to_csv(
        out_dir / 'edge_mouse_gene_expressed_in_tissue.tsv', sep='\t', index=False
    )
    logger.info(
        'Written: edge_mouse_gene_expressed_in_tissue.tsv (%d rows)', len(edge_df)
    )
    if d := edge_df.duplicated(subset=['gene_id', 'tissue_id']).sum():
        logger.warning('Duplicate gene→tissue rows in output: %d', d)

    tissue_df = (
        edge_df[['tissue_id', 'tissue_name']]
        .drop_duplicates(subset=['tissue_id'])
        .sort_values('tissue_name')
        .reset_index(drop=True)
    )
    tissue_df.to_csv(out_dir / 'node_tissue_gxd.tsv', sep='\t', index=False)
    logger.info('Written: node_tissue_gxd.tsv (%d rows)', len(tissue_df))

    logger.info('Unique genes: %d', edge_df['gene_id'].nunique())
    logger.info('Unique tissues: %d', edge_df['tissue_id'].nunique())
    logger.info(
        'TPM level distribution:\n%s',
        edge_df['tpm_level'].value_counts().to_string(),
    )

    return {'edges': edge_df, 'tissues': tissue_df}
