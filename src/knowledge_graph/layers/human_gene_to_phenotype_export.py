"""
src/knowledge_graph/layers/human_gene_to_phenotype_export.py - gene-to-phenotype export pipeline.

Builds human gene→MP top-term edges (via HPO→MP anchor mapping) and mouse
gene→MP top-term edges (via MGI PhenoGenoMP and pronto MP ontology).
Both species branches write BioCypher-ready TSVs.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import pandas as pd
from src.knowledge_graph.gene_ids import add_gene_id_column, validate_required_columns
from src.knowledge_graph.layers._mp_ontology import (
    MP_ROOT,
    collapse_to_top,
    top_children,
)
from src.knowledge_graph.layers.mouse_gene_to_phenotype_export import (
    _build_mouse_gene_mp_top_edges,
)

logger = logging.getLogger(__name__)


def load_mp_ontology(data_dir: Path):
    """
    Load the MP ontology from mp.obo in data_dir using pronto.

    Args:
        data_dir: Directory containing mp.obo.

    Returns:
        A pronto.Ontology instance.

    Raises:
        FileNotFoundError: If mp.obo is absent.
    """
    from pronto import Ontology

    mp_path = data_dir / 'mp.obo'
    if not mp_path.exists():
        raise FileNotFoundError(mp_path)
    logger.debug('Loading MP ontology from %s', mp_path)
    return Ontology(str(mp_path))


def build_top_edges(
    terms: pd.Series,
    ont,
    top_set: set[str],
    root_id: str,
    term_col: str,
    top_col: str,
) -> pd.DataFrame:
    """
    Map a series of MP term IDs to their MP top-level ancestors.

    Args:
        terms: Series of MP term ID strings.
        ont: pronto.Ontology instance.
        top_set: Set of MP top-level term IDs.
        root_id: Root term ID to exclude from results.
        term_col: Column name for the source term.
        top_col: Column name for the top-level ancestor.

    Returns:
        DataFrame with columns [term_col, top_col].
    """
    edges = []
    for tid in sorted(set(t for t in terms.dropna().unique() if isinstance(t, str))):
        for top in collapse_to_top(ont, tid, top_set, root_id):
            edges.append({term_col: tid, top_col: top})
    return pd.DataFrame(edges) if edges else pd.DataFrame(columns=[term_col, top_col])


def read_genes_to_phenotype(path: str | Path) -> pd.DataFrame:
    """
    Read the HPO genes_to_phenotype.txt file.

    Args:
        path: Path to genes_to_phenotype.txt.

    Returns:
        DataFrame with at least columns: gene_symbol, hpo_id.

    Raises:
        ValueError: If required columns are absent.
    """
    logger.debug('Reading genes_to_phenotype from %s', path)
    df = pd.read_csv(path, sep='\t', dtype=str)
    required = {'hpo_id', 'gene_symbol'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f'genes_to_phenotype missing columns: {missing}. Got: {list(df.columns)}'
        )

    n_before = len(df)
    df = df[df['gene_symbol'].notna() & (df['gene_symbol'] != '-')]
    dropped = n_before - len(df)
    if dropped:
        logger.warning(
            "Dropped %d rows with gene_symbol='-' (no known symbol in HPO)", dropped
        )
    logger.debug(
        'g2p: %d rows, %d genes, %d unique HP terms',
        len(df),
        df['gene_symbol'].nunique(),
        df['hpo_id'].nunique(),
    )
    return df


def read_hp_to_mp_top(mapping_dir: Path) -> pd.DataFrame:
    """
    Load the HP→MP top-term mapping produced by phenotype_mapping.py.

    Args:
        mapping_dir: Directory containing edge_hp_to_mp_top.tsv.

    Returns:
        DataFrame with columns: hp_id, mp_top_id, hp_to_anchor_distance.

    Raises:
        FileNotFoundError: If edge_hp_to_mp_top.tsv is absent.
        ValueError: If required columns are missing.
    """
    path = mapping_dir / 'edge_hp_to_mp_top.tsv'
    if not path.exists():
        raise FileNotFoundError(
            f'Mapping file not found: {path}\n'
            'Run phenotype_mapping.py first to generate this file.'
        )
    df = pd.read_csv(path, sep='\t', dtype=str)
    required = {'hp_id', 'mp_top_id', 'hp_to_anchor_distance'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f'edge_hp_to_mp_top.tsv missing columns: {missing}')
    logger.info(
        'HP→MP top mapping loaded: %d edges, %d unique HP terms, %d unique MP terms',
        len(df),
        df['hp_id'].nunique(),
        df['mp_top_id'].nunique(),
    )
    return df


def read_mgi_phenogenomp(path: str | Path) -> pd.DataFrame:
    """
    Read the MGI PhenoGenoMP.rpt file.

    Args:
        path: Path to MGI_PhenoGenoMP.rpt.

    Returns:
        DataFrame with columns: genotype, allele, background, mp_id, pmid,
        marker_mgi, genotype_mgi.
    """
    logger.debug('Reading MGI PhenoGenoMP from %s', path)
    cols = [
        'genotype',
        'allele',
        'background',
        'mp_id',
        'pmid',
        'marker_mgi',
        'genotype_mgi',
    ]
    df = pd.read_csv(path, sep='\t', header=None, names=cols, dtype=str)
    return df


def build_human_gene_mp_top(g2p, hp_to_mp_top):
    """
    Merge HPO gene-phenotype annotations with HP→MP top mapping.

    Args:
        g2p: genes_to_phenotype DataFrame with gene_symbol and hpo_id.
        hp_to_mp_top: HP→MP top mapping DataFrame.

    Returns:
        DataFrame with columns: gene_id, gene_symbol, mp_top_id,
        source_hp_ids, min_hp_to_anchor_dist.
    """
    validate_required_columns(
        g2p, {'gene_symbol', 'hpo_id'}, 'human phenotype export input'
    )
    validate_required_columns(
        hp_to_mp_top,
        {'hp_id', 'mp_top_id', 'hp_to_anchor_distance'},
        'HP to MP top mapping',
    )
    gene_hp = g2p[['gene_symbol', 'hpo_id']].drop_duplicates().dropna()
    merged = gene_hp.merge(
        hp_to_mp_top, left_on='hpo_id', right_on='hp_id', how='inner'
    )

    agg = (
        merged.groupby(['gene_symbol', 'mp_top_id'])
        .agg(
            source_hp_ids=('hpo_id', lambda x: '|'.join(sorted(x.unique()))),
            min_hp_to_anchor_dist=('hp_to_anchor_distance', 'min'),
        )
        .reset_index()
    )
    # MVP convention: HGNC:<symbol> is not a canonical HGNC identifier.
    # TODO: migrate human gene nodes to true HGNC IDs when available.
    agg = add_gene_id_column(agg, 'human')
    agg = agg[
        [
            'gene_id',
            'gene_symbol',
            'mp_top_id',
            'source_hp_ids',
            'min_hp_to_anchor_dist',
        ]
    ]
    logger.debug('Human gene→MP top edges: %d', len(agg))
    return agg


def run_pipeline(
    genes_to_phenotype_path: str | Path,
    mgi_path: str | Path,
    data_dir: str | Path,
    mapping_dir: str | Path,
    gene_info_path: str | Path,
    out_dir: str | Path = './out',
    species: str = 'both',
) -> dict[str, pd.DataFrame]:
    """
    Run the gene-to-phenotype export pipeline for human and/or mouse.

    Args:
        genes_to_phenotype_path: Path to HPO genes_to_phenotype.txt.
        mgi_path: Path to MGI_PhenoGenoMP.rpt.
        data_dir: Directory containing mp.obo (and optionally hp.obo).
        mapping_dir: Directory with edge_hp_to_mp_top.tsv from phenotype_mapping.py.
        gene_info_path: Path to NCBI gene_info.gz.
        out_dir: Output directory for TSV files.
        species: 'human', 'mouse', or 'both'.

    Returns:
        Dict with keys: human_gene_mp_top, mouse_gene_mp, mp_in_top_from_mgi,
        mouse_gene_mp_top.

    Raises:
        RuntimeError: If either species branch fails.
    """
    out_dir = Path(out_dir)
    data_dir = Path(data_dir)
    mapping_dir = Path(mapping_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info('Gene-to-phenotype export pipeline')
    logger.info('genes_to_phenotype: %s', genes_to_phenotype_path)
    logger.info('mapping_dir: %s', mapping_dir)
    logger.info('data_dir: %s', data_dir)
    logger.info('out_dir: %s', out_dir)
    logger.info('species: %s', species)

    build_human = species in ('human', 'both')
    build_mouse = species in ('mouse', 'both')

    names_src = mapping_dir / 'node_mp_top_names.tsv'
    names_dst = out_dir / 'node_mp_top_names.tsv'
    # Use resolve() to handle symlinks correctly before comparing paths (§13 fix).
    if names_src.exists() and names_src.resolve() != names_dst.resolve():
        shutil.copy(names_src, names_dst)
        logger.debug('Copied node_mp_top_names.tsv to out_dir')

    human_gene_mp_top = pd.DataFrame()
    if build_human:
        logger.info('--- Human branch ---')
        try:
            g2p = read_genes_to_phenotype(genes_to_phenotype_path)
            hp_to_mp_top = read_hp_to_mp_top(mapping_dir)

            human_gene_mp_top = build_human_gene_mp_top(g2p, hp_to_mp_top)
            human_gene_mp_top.to_csv(
                out_dir / 'edge_human_gene_has_mp_top.tsv', sep='\t', index=False
            )
            logger.info(
                'Written: edge_human_gene_has_mp_top.tsv (%d rows)',
                len(human_gene_mp_top),
            )
            _dupes = human_gene_mp_top.duplicated(subset=['gene_id', 'mp_top_id']).sum()
            if _dupes:
                logger.warning('Duplicate gene→MP top rows in output: %d', _dupes)

            mapped_hp = set(hp_to_mp_top['hp_id'].dropna())
            unmapped = sorted(set(g2p['hpo_id'].dropna()) - mapped_hp)
            pd.DataFrame({'unmapped_hp_id': unmapped}).to_csv(
                out_dir / 'qc_unmapped_hp.tsv', sep='\t', index=False
            )
            logger.info('Human gene→MP top edges: %d', len(human_gene_mp_top))
            logger.info(
                'Unmapped HP terms: %d / %d',
                len(unmapped),
                g2p['hpo_id'].nunique(),
            )
        except Exception as exc:
            raise RuntimeError(f'Human branch failed: {exc}') from exc

    mouse_gene_mp = mp_in_top_from_mgi = mouse_gene_mp_top = pd.DataFrame()
    if build_mouse:
        logger.info('--- Mouse branch ---')
        try:
            from src.utils.gene_info import build_mouse_symbol_maps

            _, mgi_to_symbol = build_mouse_symbol_maps(gene_info_path)

            mp_ont = load_mp_ontology(data_dir)
            mp_top = top_children(mp_ont, MP_ROOT)

            mp_top_names = pd.DataFrame(
                [
                    {
                        'mp_id': mp_id,
                        'name': (mp_ont[mp_id].name or '') if mp_id in mp_ont else '',
                    }
                    for mp_id in sorted(mp_top)
                ]
            )
            mp_top_names.to_csv(
                out_dir / 'node_mp_top_names.tsv', sep='\t', index=False
            )
            logger.debug('Written: node_mp_top_names.tsv (%d terms)', len(mp_top_names))

            mgi = read_mgi_phenogenomp(mgi_path)

            raw_mouse = (
                mgi[['marker_mgi', 'mp_id']]
                .dropna(subset=['marker_mgi', 'mp_id'])
                .drop_duplicates()
                .copy()
            )
            raw_mouse['mgi_single'] = raw_mouse['marker_mgi'].str.split('|')
            raw_mouse = raw_mouse.explode('mgi_single')
            raw_mouse['mgi_single'] = raw_mouse['mgi_single'].str.strip()
            raw_mouse['gene_symbol'] = raw_mouse['mgi_single'].map(mgi_to_symbol)

            n_total = raw_mouse['mgi_single'].nunique()
            n_mapped = raw_mouse.dropna(subset=['gene_symbol'])['mgi_single'].nunique()
            if n_total - n_mapped:
                logger.warning(
                    'QC: %d / %d MGI IDs could not be mapped to a gene symbol',
                    n_total - n_mapped,
                    n_total,
                )

            mouse_gene_mp = raw_mouse.dropna(subset=['gene_symbol'])[
                ['gene_symbol', 'mp_id']
            ].drop_duplicates()
            mouse_gene_mp = add_gene_id_column(mouse_gene_mp, 'mouse')
            mouse_gene_mp = mouse_gene_mp[['gene_id', 'gene_symbol', 'mp_id']]
            mouse_gene_mp.to_csv(
                out_dir / 'edge_mouse_gene_has_mp.tsv', sep='\t', index=False
            )

            mp_in_top_from_mgi = build_top_edges(
                mgi['mp_id'], mp_ont, mp_top, MP_ROOT, 'mp_id', 'mp_top_id'
            )
            mp_in_top_from_mgi.to_csv(
                out_dir / 'edge_mp_in_top_from_mgi.tsv', sep='\t', index=False
            )

            mouse_gene_mp_top = _build_mouse_gene_mp_top_edges(
                mouse_gene_mp, mp_in_top_from_mgi
            )
            mouse_gene_mp_top.to_csv(
                out_dir / 'edge_mouse_gene_has_mp_top.tsv', sep='\t', index=False
            )
            _dupes = mouse_gene_mp_top.duplicated(subset=['gene_id', 'mp_top_id']).sum()
            if _dupes:
                logger.warning('Duplicate mouse gene→MP top rows in output: %d', _dupes)

            logger.info('Mouse gene→MP edges: %d', len(mouse_gene_mp))
            logger.info('Mouse gene→MP top: %d', len(mouse_gene_mp_top))
        except Exception as exc:
            raise RuntimeError(f'Mouse branch failed: {exc}') from exc

    return {
        'human_gene_mp_top': human_gene_mp_top,
        'mouse_gene_mp': mouse_gene_mp,
        'mp_in_top_from_mgi': mp_in_top_from_mgi,
        'mouse_gene_mp_top': mouse_gene_mp_top,
    }


if __name__ == '__main__':
    import argparse
    import logging

    from src.utils.logging_config import setup_logging

    p = argparse.ArgumentParser(
        description=(
            'Build gene→MP system-level edges. '
            'Requires edge_hp_to_mp_top.tsv from phenotype_mapping.py.'
        )
    )
    p.add_argument('--genes-to-phenotype', required=True)
    p.add_argument('--mgi', required=True)
    p.add_argument('--data-dir', required=True, help='Directory with mp.obo')
    p.add_argument(
        '--mapping-dir',
        required=True,
        help='Directory with outputs of phenotype_mapping.py',
    )
    p.add_argument('--gene-info', required=True)
    p.add_argument('--out', default='./out')
    p.add_argument('--species', default='both', choices=['human', 'mouse', 'both'])
    p.add_argument(
        '--log-level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR']
    )
    args = p.parse_args()

    setup_logging(out_dir=args.out, console_level=getattr(logging, args.log_level))

    run_pipeline(
        genes_to_phenotype_path=args.genes_to_phenotype,
        mgi_path=args.mgi,
        data_dir=args.data_dir,
        mapping_dir=args.mapping_dir,
        gene_info_path=args.gene_info,
        out_dir=args.out,
        species=args.species,
    )
