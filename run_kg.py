"""
run_kg.py - CLI entry point for knowledge graph construction.

Runs all data-layer export pipelines (phenotype, GO, expression, PPI, Uberon)
for the requested species, then calls BioCypher to assemble the heterogeneous
graph and write Neo4j-ready import files.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from datetime import datetime
from pathlib import Path

from src.knowledge_graph.build_graph import DEFAULT_SCHEMA, build_species
from src.knowledge_graph.layers.gene_ontology_export import (
    run_go_pipeline,
    run_mouse_go_pipeline,
)
from src.knowledge_graph.layers.human_gene_expression_export import (
    run_expression_pipeline,
)
from src.knowledge_graph.layers.human_gene_to_phenotype_export import (
    run_pipeline as run_phenotype_pipeline,
)
from src.knowledge_graph.layers.mouse_gene_expression_export import (
    run_mouse_expression_pipeline,
)
from src.knowledge_graph.layers.ppi_export import HUMAN, MOUSE, run_ppi_pipeline
from src.knowledge_graph.layers.uberon_export import run_uberon_pipeline
from src.knowledge_graph.layers.uberon_mapping import build_uberon_mapping
from src.knowledge_graph.preprocessing.phenotype_mapping import run_mapping
from src.utils.logging_config import setup_logging

_PROJECT_ROOT = Path(__file__).resolve().parent
_RAW_DATA_DIR = _PROJECT_ROOT / 'data' / 'raw'
_PROCESSED_DIR = _PROJECT_ROOT / 'data' / 'processed'
_BIOCYPHER_DIR = _PROJECT_ROOT / 'data' / 'biocypher'

# ---------------------------------------------------------------------------
# Raw-file resolution
# ---------------------------------------------------------------------------

_RAW_FILE_ARGS: dict[str, str] = {
    'genes_to_phenotype': 'genes_to_phenotype.txt',
    'mgi': 'MGI_PhenoGenoMP.rpt',
    'hp2mp': 'HP2MP.tsv',
    'gene2go': 'gene2go.gz',
    'gene_info': 'gene_info.gz',
    'gtex': 'GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz',
    'gxd_tar': 'gxdrnaseq.tar.gz',
    'gxd_uberon_mapping': 'gxd_uberon_mapping.tsv',
    'biogrid': 'BIOGRID-ALL-5.0.256.tab3.txt',
    'uberon_obo': 'basic.obo',
}


def _resolve_raw_inputs(args: argparse.Namespace) -> argparse.Namespace:
    """
    Resolve raw input file paths from CLI args and default data directory.

    Args:
        args: Parsed argument namespace.

    Returns:
        Updated namespace with all file paths resolved to Path objects.
    """
    data_dir = Path(args.data_dir) if args.data_dir else _RAW_DATA_DIR
    args.data_dir = data_dir
    for arg_name, filename in _RAW_FILE_ARGS.items():
        value = getattr(args, arg_name)
        path = Path(value) if value else data_dir / filename
        setattr(args, arg_name, path)
    return args


# ---------------------------------------------------------------------------
# Data-file preflight check
# ---------------------------------------------------------------------------

_REQUIRED_FILES: list[tuple[str, str, str]] = [
    ('hp.obo', 'Human Phenotype Ontology', 'both'),
    ('mp.obo', 'Mammalian Phenotype Ontology', 'both'),
    ('genes_to_phenotype.txt', 'HPO gene-to-phenotype annotations', 'human'),
    ('HP2MP.tsv', 'HP→MP anchor mapping from MGI', 'human'),
    ('MGI_PhenoGenoMP.rpt', 'MGD genotype-to-MP annotations', 'mouse'),
    ('gene_info.gz', 'NCBI gene info (symbol + MGI xrefs)', 'mouse'),
    ('gene2go.gz', 'NCBI Gene-to-GO mappings', 'both'),
    (
        'GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz',
        'GTEx median gene expression (GCT)',
        'human',
    ),
    ('gxdrnaseq.tar.gz', 'MGI GXD RNA-seq archive', 'mouse'),
    ('BIOGRID-ALL-5.0.256.tab3.txt', 'BioGRID protein-protein interactions', 'both'),
    ('basic.obo', 'Uberon anatomical ontology (basic)', 'human'),
]


def _file_needed(species: str, req_species: str) -> bool:
    """
    Return True if a file with the given req_species is needed for the run species.

    Args:
        species: Target species for the current run ('human', 'mouse', or 'both').
        req_species: Species requirement for the file ('human', 'mouse', or 'both').

    Returns:
        False only when the species and req_species are opposite single-species values.
    """
    if req_species == 'human' and species == 'mouse':
        return False
    if req_species == 'mouse' and species == 'human':
        return False
    return True


def _human_readable_size(nbytes: int) -> str:
    """
    Format a byte count as a human-readable string (e.g. '1.4M').

    Args:
        nbytes: File size in bytes.

    Returns:
        Human-readable size string with unit suffix.
    """
    for unit in ('B', 'K', 'M', 'G', 'T'):
        if abs(nbytes) < 1024 or unit == 'T':
            return f'{nbytes:.1f}{unit}' if unit != 'B' else f'{nbytes}B'
        nbytes //= 1024
    return f'{nbytes:.1f}T'


def check_data(args: argparse.Namespace) -> None:
    """
    Verify all required input files are present; exit with an error if any are missing.

    Args:
        args: Resolved argument namespace (after _resolve_raw_inputs).
    """
    data_dir = Path(args.data_dir)
    path_overrides: dict[str, Path] = {
        'genes_to_phenotype.txt': Path(args.genes_to_phenotype),
        'MGI_PhenoGenoMP.rpt': Path(args.mgi),
        'HP2MP.tsv': Path(args.hp2mp),
        'gene2go.gz': Path(args.gene2go),
        'gene_info.gz': Path(args.gene_info),
        'GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz': Path(args.gtex),
        'gxdrnaseq.tar.gz': Path(args.gxd_tar),
        'BIOGRID-ALL-5.0.256.tab3.txt': Path(args.biogrid),
        'basic.obo': Path(args.uberon_obo)
    }

    print(f'Checking required data files in: {data_dir}/')
    print('-' * 70)

    missing: list[tuple[str, str]] = []
    checked = 0
    for filename, description, req_species in _REQUIRED_FILES:
        if not _file_needed(args.species, req_species):
            continue
        filepath = path_overrides.get(filename, data_dir / filename)
        checked += 1
        if filepath.exists():
            size = _human_readable_size(filepath.stat().st_size)
            print(f'  [OK]  {filename}  ({size})  {description}')
        else:
            print(f'  [MISSING]  {filename}  {description}')
            missing.append((str(filepath), description))

    print('-' * 70)
    if not missing:
        print(f'All {checked} files present. Ready to run the pipeline.\n')
        return

    print()
    print('=' * 70)
    print('MISSING DATA FILES')
    print('=' * 70)
    for filepath, description in missing:
        print(f'  {filepath}\n    {description}')
    print(
        '\nDownload missing raw files into data/raw:\n'
        '  uv run python scripts/download_data.py\n\n'
        'Then re-run this command.'
    )
    print('=' * 70)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Per-species export
# ---------------------------------------------------------------------------


def _export_human(args: argparse.Namespace, out_dir: Path) -> None:
    """
    Run all human data-layer export pipelines.

    Args:
        args: Resolved CLI argument namespace.
        out_dir: Output directory for processed TSV files.
    """
    print('=' * 70)
    print('Step 0 - HP→MP system-level mapping')
    print('=' * 70)
    run_mapping(
        genes_to_phenotype_path=args.genes_to_phenotype,
        hp2mp_path=args.hp2mp,
        data_dir=args.data_dir,
        out_dir=out_dir,
        exclude_mortality=not args.include_mortality,
    )

    print('=' * 70)
    print('Step 1 - Phenotype edges (HPO)')
    print('=' * 70)
    run_phenotype_pipeline(
        genes_to_phenotype_path=args.genes_to_phenotype,
        mgi_path=args.mgi,
        data_dir=args.data_dir,
        mapping_dir=out_dir,
        gene_info_path=args.gene_info,
        out_dir=out_dir,
        species='human',
    )

    print('=' * 70)
    print('Step 2 - GO annotation edges (taxon 9606)')
    print('=' * 70)
    run_go_pipeline(
        gene2go_path=args.gene2go,
        genes_to_phenotype_path=args.genes_to_phenotype,
        out_dir=out_dir,
    )

    print('=' * 70)
    print('Step 3 - GTEx expression edges')
    print('=' * 70)
    run_expression_pipeline(
        gtex_path=args.gtex,
        out_dir=out_dir,
        tpm_threshold=args.tpm_threshold,
    )

    print('=' * 70)
    print('Step 4 - BioGRID protein-protein interactions (human)')
    print('=' * 70)
    run_ppi_pipeline(
        biogrid_path=args.biogrid,
        cfg=HUMAN,
        out_dir=out_dir,
        min_publications=args.min_publications,
    )

    print('=' * 70)
    print('Step 5 - Uberon anatomical ontology')
    print('=' * 70)
    expression_tsv = out_dir / 'edge_human_gene_expressed_in_tissue.tsv'
    run_uberon_pipeline(
        out_dir=out_dir,
        uberon_obo_path=args.uberon_obo,
        expression_tsv=expression_tsv if expression_tsv.exists() else None,
    )


def _export_mouse(args: argparse.Namespace, out_dir: Path) -> None:
    """
    Run all mouse data-layer export pipelines.

    Args:
        args: Resolved CLI argument namespace.
        out_dir: Output directory for processed TSV files.
    """
    print('=' * 70)
    print('Step 1 - Phenotype edges (MGD)')
    print('=' * 70)
    run_phenotype_pipeline(
        genes_to_phenotype_path=args.genes_to_phenotype,
        mgi_path=args.mgi,
        data_dir=args.data_dir,
        mapping_dir=out_dir,
        gene_info_path=args.gene_info,
        out_dir=out_dir,
        species='mouse',
    )

    print('=' * 70)
    print('Step 2 - GO annotation edges (taxon 10090)')
    print('=' * 70)
    run_mouse_go_pipeline(
        gene2go_path=args.gene2go,
        gene_info_path=args.gene_info,
        out_dir=out_dir,
    )

    print('=' * 70)
    print('Step 3 - GXD expression edges')
    print('=' * 70)
    if not args.gxd_uberon_mapping.exists():
        print('gxd_uberon_mapping.tsv not found - generating from basic.obo...')
        build_uberon_mapping(
            tar_path=args.gxd_tar,
            obo_path=args.uberon_obo,
            out_path=args.gxd_uberon_mapping,
        )
        print(
            f'Generated: {args.gxd_uberon_mapping}\n'
            'Review unmatched rows (empty uberon_id) and re-run if needed.'
        )
    run_mouse_expression_pipeline(
        tar_path=args.gxd_tar,
        data_dir=args.data_dir,
        out_dir=out_dir,
        uberon_mapping_path=args.gxd_uberon_mapping,
    )

    print('=' * 70)
    print('Step 4 - BioGRID protein-protein interactions (mouse)')
    print('=' * 70)
    run_ppi_pipeline(
        biogrid_path=args.biogrid,
        cfg=MOUSE,
        out_dir=out_dir,
        min_publications=args.min_publications,
    )


_EXPORT_FN = {
    'human': _export_human,
    'mouse': _export_mouse,
}

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    """
    Parse and resolve CLI arguments for the knowledge graph build pipeline.

    Returns:
        Resolved argparse.Namespace with all file paths as Path objects.
    """
    p = argparse.ArgumentParser(
        description='Build gene-phenotype BioCypher knowledge graphs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Raw input overrides (all default to <data-dir>/<standard-filename>)
    p.add_argument(
        '--data-dir',
        default=None,
        help=f'Directory containing raw ontology/source files (default: {_RAW_DATA_DIR})',
    )
    p.add_argument('--genes-to-phenotype', default=None)
    p.add_argument('--mgi', default=None)
    p.add_argument('--hp2mp', default=None)
    p.add_argument('--gene2go', default=None)
    p.add_argument('--gene-info', default=None)
    p.add_argument('--gtex', default=None)
    p.add_argument('--gxd-tar', default=None)
    p.add_argument('--gxd-uberon-mapping', default=None)
    p.add_argument('--biogrid', default=None)
    p.add_argument('--uberon-obo', default=None)

    # Filtering
    p.add_argument(
        '--tpm-threshold',
        type=float,
        default=5.0,
        help='Minimum median TPM for a gene-tissue edge',
    )
    p.add_argument(
        '--min-publications',
        type=int,
        default=1,
        metavar='N',
        help='Minimum distinct PubMed IDs per PPI pair',
    )

    # Species and phenotype options
    p.add_argument('--species', choices=['human', 'mouse', 'both'], default='both')
    p.add_argument(
        '--include-mortality',
        action='store_true',
        help='Keep MP:0010768 (mortality/aging) in HP→MP mapping',
    )

    # BioCypher schema overrides (rarely needed)
    p.add_argument(
        '--schema-config',
        default=None,
        help='Explicit schema YAML (single-species builds only)',
    )
    p.add_argument('--biocypher-config', default=None)

    # Misc
    p.add_argument(
        '--skip-export',
        action='store_true',
        help='Skip TSV export; use existing files in data/processed/',
    )

    return _resolve_raw_inputs(p.parse_args())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """
    Orchestrate the full knowledge graph build: export TSVs then run BioCypher.
    """
    args = parse_args()
    setup_logging(out_dir=_PROJECT_ROOT / 'logs')

    targets = ['human', 'mouse'] if args.species == 'both' else [args.species]
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    if not args.skip_export:
        check_data(args)

    for species in targets:
        processed_dir = _PROCESSED_DIR / species
        biocypher_dir = _BIOCYPHER_DIR / f'{species}_{timestamp}'

        print(f'\n{"=" * 70}')
        print(f'  {species.upper()}')
        print(f'  processed → {processed_dir}')
        print(f'  biocypher → {biocypher_dir}')
        print(f'{"=" * 70}\n')

        if not args.skip_export:
            # Always start with a clean processed dir
            if processed_dir.exists():
                shutil.rmtree(processed_dir)
            processed_dir.mkdir(parents=True)

            _EXPORT_FN[species](args, processed_dir)
            print(f'\n→ TSVs written to {processed_dir.resolve()}\n')

        schema_path = (
            Path(args.schema_config) if args.schema_config else DEFAULT_SCHEMA[species]
        )
        bc_config = Path(args.biocypher_config) if args.biocypher_config else None

        build_species(
            species=species,
            data_dir=processed_dir,
            out_dir=biocypher_dir,
            schema_config_path=schema_path,
            biocypher_config_path=bc_config,
            skip_layers=[],
            allow_missing_layers=False,
        )
        print(f'\n→ Graph written to {biocypher_dir.resolve()}\n')


if __name__ == '__main__':
    main()
