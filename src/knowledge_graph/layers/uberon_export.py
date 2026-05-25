"""
src/knowledge_graph/layers/uberon_export.py - Uberon anatomical ontology export pipeline.

Maps curated GTEx tissue IDs to Uberon terms, optionally populates term
metadata from an OBO file, and writes node/edge TSVs for the human Uberon layer.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
from src.knowledge_graph.layers.human_gene_expression_export import _gtex_tissue_id

logger = logging.getLogger(__name__)


GTEX_TO_UBERON: dict[str, str] = {
    'GTEX:Brain_Amygdala': 'UBERON:0001876',
    'GTEX:Brain_Anterior_cingulate_cortex_BA24': 'UBERON:0006101',
    'GTEX:Brain_Caudate_basal_ganglia': 'UBERON:0001873',
    'GTEX:Brain_Cerebellar_Hemisphere': 'UBERON:0002037',
    'GTEX:Brain_Cerebellum': 'UBERON:0002037',
    'GTEX:Brain_Cortex': 'UBERON:0001851',
    'GTEX:Brain_Frontal_Cortex_BA9': 'UBERON:0001870',
    'GTEX:Brain_Hippocampus': 'UBERON:0002310',
    'GTEX:Brain_Hypothalamus': 'UBERON:0001898',
    'GTEX:Brain_Nucleus_accumbens_basal_ganglia': 'UBERON:0001882',
    'GTEX:Brain_Putamen_basal_ganglia': 'UBERON:0001874',
    'GTEX:Brain_Spinal_cord_cervical_c_1': 'UBERON:0002726',
    'GTEX:Brain_Substantia_nigra': 'UBERON:0002038',
    'GTEX:Heart_Atrial_Appendage': 'UBERON:0006618',
    'GTEX:Heart_Left_Ventricle': 'UBERON:0002084',
    'GTEX:Artery_Aorta': 'UBERON:0000947',
    'GTEX:Artery_Coronary': 'UBERON:0001621',
    'GTEX:Artery_Tibial': 'UBERON:0007610',
    'GTEX:Colon_Sigmoid': 'UBERON:0001159',
    'GTEX:Colon_Transverse': 'UBERON:0001157',
    'GTEX:Small_Intestine_Terminal_Ileum': 'UBERON:0002116',
    'GTEX:Stomach': 'UBERON:0000945',
    'GTEX:Liver': 'UBERON:0002107',
    'GTEX:Pancreas': 'UBERON:0001264',
    'GTEX:Esophagus_Gastroesophageal_Junction': 'UBERON:0007650',
    'GTEX:Esophagus_Mucosa': 'UBERON:0002469',
    'GTEX:Esophagus_Muscularis': 'UBERON:0004648',
    'GTEX:Ovary': 'UBERON:0000992',
    'GTEX:Uterus': 'UBERON:0000995',
    'GTEX:Vagina': 'UBERON:0000996',
    'GTEX:Prostate': 'UBERON:0002367',
    'GTEX:Testis': 'UBERON:0000473',
    'GTEX:Adrenal_Gland': 'UBERON:0002369',
    'GTEX:Thyroid': 'UBERON:0002046',
    'GTEX:Pituitary': 'UBERON:0000007',
    'GTEX:Lung': 'UBERON:0002048',
    'GTEX:Muscle_Skeletal': 'UBERON:0001134',
    'GTEX:Minor_Salivary_Gland': 'UBERON:0006330',
    'GTEX:Skin_Not_Sun_Exposed_Suprapubic': 'UBERON:0000014',
    'GTEX:Skin_Sun_Exposed_Lower_leg': 'UBERON:0000014',
    'GTEX:Nerve_Tibial': 'UBERON:0001323',
    'GTEX:Whole_Blood': 'UBERON:0000178',
    'GTEX:Spleen': 'UBERON:0002106',
    'GTEX:Kidney_Cortex': 'UBERON:0001225',
    'GTEX:Kidney_Medulla': 'UBERON:0001293',
    'GTEX:Bladder': 'UBERON:0001255',
    'GTEX:Adipose_Subcutaneous': 'UBERON:0002190',
    'GTEX:Adipose_Visceral_Omentum': 'UBERON:0003688',
    'GTEX:Breast_Mammary_Tissue': 'UBERON:0001911',
    'GTEX:Cells_Cultured_fibroblasts': 'UBERON:0002082',
    'GTEX:Cells_EBV_transformed_lymphocytes': 'UBERON:0000029',
    'GTEX:Cervix_Ectocervix': 'UBERON:0012249',
    'GTEX:Cervix_Endocervix': 'UBERON:0000458',
    'GTEX:Fallopian_Tube': 'UBERON:0003889',
}


def _parse_uberon_obo(obo_path: str | Path) -> pd.DataFrame:
    """
    Parse a Uberon OBO file and return a DataFrame of UBERON term metadata.

    Args:
        obo_path: Path to the Uberon basic OBO file.

    Returns:
        DataFrame with columns: uberon_id, name, definition, synonyms.

    Raises:
        ImportError: If obonet is not installed.
    """
    try:
        import obonet
    except ImportError as exc:
        raise ImportError(
            'obonet is required to parse the Uberon OBO file.\n'
            'Install it with:  uv add obonet'
        ) from exc

    logger.info('Parsing Uberon OBO: %s', obo_path)
    graph = obonet.read_obo(str(obo_path))
    logger.info('OBO graph loaded - %d total terms', graph.number_of_nodes())

    rows: list[dict] = []
    for node_id, data in graph.nodes(data=True):
        if not node_id.startswith('UBERON:'):
            continue
        name = data.get('name', '')

        raw_def: str = data.get('def', '')
        definition = raw_def.split('" [')[0].lstrip('"') if raw_def else ''

        raw_syns: list[str] = data.get('synonym', [])
        synonym_names = [s.split('"')[1] for s in raw_syns if '"' in s]
        rows.append(
            {
                'uberon_id': node_id,
                'name': name,
                'definition': definition,
                'synonyms': '|'.join(synonym_names),
            }
        )

    df = pd.DataFrame(rows, columns=['uberon_id', 'name', 'definition', 'synonyms'])
    logger.info('Uberon terms parsed: %d', len(df))
    return df


def _parse_uberon_obo_fallback(uberon_ids: set[str]) -> pd.DataFrame:
    """
    Build stub Uberon term nodes with empty metadata when no OBO file is available.

    Args:
        uberon_ids: Set of UBERON IDs referenced by the tissue mapping.

    Returns:
        DataFrame with columns: uberon_id, name, definition, synonyms (all empty).
    """
    logger.warning(
        'No Uberon OBO file provided - writing stub nodes with empty metadata. '
        'Re-run with --uberon-obo <path> to populate names and definitions.'
    )
    rows = [
        {'uberon_id': uid, 'name': '', 'definition': '', 'synonyms': ''}
        for uid in sorted(uberon_ids)
    ]
    return pd.DataFrame(rows, columns=['uberon_id', 'name', 'definition', 'synonyms'])


def build_tissue_to_uberon_edges(
    mapping: dict[str, str],
    tissue_ids_in_graph: set[str] | None = None,
) -> pd.DataFrame:
    """
    Build tissue→Uberon edge rows from the curated GTEx→Uberon mapping.

    Args:
        mapping: Dict mapping GTEX tissue ID → UBERON ID.
        tissue_ids_in_graph: Optional set of tissue IDs to restrict output to.

    Returns:
        DataFrame with columns: gtex_tissue_id, uberon_id, gtex_tissue_name.
    """
    rows: list[dict] = []
    for gtex_id, uberon_id in mapping.items():
        if tissue_ids_in_graph is not None and gtex_id not in tissue_ids_in_graph:
            logger.debug('Tissue not in expression layer, skipping: %s', gtex_id)
            continue
        tissue_name = gtex_id.replace('GTEX:', '').replace('_', ' ')
        rows.append(
            {
                'gtex_tissue_id': gtex_id,
                'uberon_id': uberon_id,
                'gtex_tissue_name': tissue_name,
            }
        )

    df = pd.DataFrame(rows, columns=['gtex_tissue_id', 'uberon_id', 'gtex_tissue_name'])
    logger.info('Tissue→Uberon edges built: %d', len(df))
    return df


def _read_expression_tissue_ids(expression_tsv: str | Path) -> set[str]:
    """
    Extract tissue IDs from an expression edge TSV for filtering.

    Args:
        expression_tsv: Path to edge_human_gene_expressed_in_tissue.tsv.

    Returns:
        Set of tissue ID strings, or empty set on missing file or columns.
    """
    path = Path(expression_tsv)
    if not path.exists():
        logger.warning(
            'Expression TSV not found: %s - using full curated mapping', path
        )
        return set()

    df = pd.read_csv(path, sep='\t', dtype=str)
    if 'tissue_id' in df.columns:
        ids = set(df['tissue_id'].dropna().unique())
    elif 'tissue_name' in df.columns:
        ids = {_gtex_tissue_id(n) for n in df['tissue_name'].dropna().unique()}
    else:
        logger.warning(
            "Expression TSV has neither 'tissue_id' nor 'tissue_name' column - "
            'using full curated mapping.'
        )
        return set()

    logger.debug('Tissue IDs found in expression layer: %d', len(ids))
    return ids


def run_uberon_pipeline(
    out_dir: str | Path,
    uberon_obo_path: str | Path | None = None,
    expression_tsv: str | Path | None = None,
    mapping: dict[str, str] | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Run the Uberon anatomical ontology export pipeline.

    Args:
        out_dir: Output directory for TSV files.
        uberon_obo_path: Optional path to the Uberon basic OBO file.
        expression_tsv: Optional path to expression edge TSV for tissue filtering.
        mapping: GTEx→Uberon mapping dict; defaults to GTEX_TO_UBERON if None.

    Returns:
        Dict with keys 'uberon_terms' and 'tissue_uberon_edges'.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mapping = mapping or GTEX_TO_UBERON

    logger.info('=' * 70)
    logger.info('Starting Uberon pipeline')
    logger.info('Output dir: %s', out_dir)
    logger.info('OBO file: %s', uberon_obo_path or '<not provided - stub mode>')
    logger.info('Expression TSV: %s', expression_tsv or '<not provided>')
    logger.info('Mapping size: %d GTEx tissues', len(mapping))
    logger.info('=' * 70)

    tissue_ids_in_graph: set[str] | None = None
    if expression_tsv:
        tissue_ids_in_graph = _read_expression_tissue_ids(expression_tsv)
        if not tissue_ids_in_graph:
            logger.warning('Expression TSV yielded 0 tissue IDs - using full mapping.')
            tissue_ids_in_graph = None

    tissue_uberon = build_tissue_to_uberon_edges(mapping, tissue_ids_in_graph)
    tissue_uberon_path = out_dir / 'edge_tissue_mapped_to_uberon.tsv'
    tissue_uberon.to_csv(tissue_uberon_path, sep='\t', index=False)
    logger.info('Written: %s  (%d rows)', tissue_uberon_path.name, len(tissue_uberon))

    referenced_uberon_ids = set(tissue_uberon['uberon_id'].dropna().unique())
    logger.debug('Unique Uberon IDs referenced: %d', len(referenced_uberon_ids))

    if uberon_obo_path and Path(uberon_obo_path).exists():
        all_uberon = _parse_uberon_obo(uberon_obo_path)
        uberon_terms = all_uberon[
            all_uberon['uberon_id'].isin(referenced_uberon_ids)
        ].copy()
        missing = referenced_uberon_ids - set(uberon_terms['uberon_id'])
        if missing:
            logger.warning(
                '%d referenced Uberon IDs not found in the OBO file: %s',
                len(missing),
                sorted(missing)[:10],
            )
        logger.info(
            'Uberon term nodes: %d / %d referenced IDs found in OBO',
            len(uberon_terms),
            len(referenced_uberon_ids),
        )
    else:
        uberon_terms = _parse_uberon_obo_fallback(referenced_uberon_ids)

    uberon_node_path = out_dir / 'node_uberon_terms.tsv'
    uberon_terms.to_csv(uberon_node_path, sep='\t', index=False)
    logger.info('Written: %s  (%d rows)', uberon_node_path.name, len(uberon_terms))

    all_mapped = set(tissue_uberon['gtex_tissue_id'])
    all_curated = set(mapping.keys())
    unmapped = sorted(all_curated - all_mapped)
    if unmapped:
        qc_path = out_dir / 'qc_uberon_unmapped_tissues.tsv'
        pd.DataFrame({'gtex_tissue_id': unmapped}).to_csv(
            qc_path, sep='\t', index=False
        )
        logger.info(
            'QC: %d curated GTEx tissues were not written (filtered out or '
            'missing in expression layer) → %s',
            len(unmapped),
            qc_path.name,
        )

    logger.info('Uberon pipeline complete.')
    return {'uberon_terms': uberon_terms, 'tissue_uberon_edges': tissue_uberon}


if __name__ == '__main__':
    from src.utils.logging_config import setup_logging

    p = argparse.ArgumentParser(
        description='Export Uberon anatomical ontology edges for the gene-phenotype KG.'
    )
    p.add_argument(
        '--uberon-obo',
        default=None,
        help=(
            'Path to the Uberon basic OBO file '
            '(uberon/basic.obo, download from http://purl.obolibrary.org/obo/uberon/basic.obo). '
            'When omitted, stub nodes with empty metadata are written.'
        ),
    )
    p.add_argument(
        '--expression-tsv',
        default=None,
        help=(
            'Path to edge_human_gene_expressed_in_tissue.tsv (or equivalent). '
            'Used to restrict the output to tissues present in the expression layer. '
            'When omitted, all curated GTEx→Uberon mappings are written.'
        ),
    )
    p.add_argument(
        '--out',
        default='./out',
        help='Output directory for TSV files (default: ./out).',
    )
    args = p.parse_args()

    setup_logging(out_dir=args.out)

    run_uberon_pipeline(
        out_dir=args.out,
        uberon_obo_path=args.uberon_obo,
        expression_tsv=args.expression_tsv,
    )
