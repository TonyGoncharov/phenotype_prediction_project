"""
run.py - single entry point for the gene-phenotype knowledge graph pipeline.

Usage
-----
  uv run python run.py                          # full pipeline, both graphs
  uv run python run.py --species human
  uv run python run.py --skip-export            # skip TSV generation
  uv run python run.py --clean-tsv --clean-graph
  uv run python run.py --help

  uv run python run.py --species human --clean-tsv --clean-graph

Required data files (put in data/)
-------------------------------------
  HP2MP.tsv, genes_to_phenotype.txt, MGI_PhenoGenoMP.rpt, gene2go.gz,
  gene_info.gz, hp.obo, mp.obo, BIOGRID-ALL-*.tab3.txt, uberon_basic.obo
  (see --help for per-file paths and sources)
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from src.layers.phenotype_mapping import run_mapping
from src.layers.gene_to_phenotype_export import run_pipeline as run_phenotype_pipeline
from src.layers.gene_ontology_export import run_go_pipeline, run_mouse_go_pipeline
from src.layers.gene_expression_export import run_expression_pipeline
from src.layers.ppi_export import run_ppi_pipeline
from src.layers.uberon_export import run_uberon_pipeline
from src.build_graph import build
from src.utils.logging_config import setup_logging


# ── Data-file validation ──────────────────────────────────────────── #

# (filename, description, needed_for_species, layer)
# species: "human", "mouse", or "both"
# layer:   None = always needed; "go" = only when GO layer is active
_REQUIRED_FILES: list[tuple[str, str, str, str | None]] = [
    ("hp.obo",                    "Human Phenotype Ontology",             "both",  None),
    ("mp.obo",                    "Mammalian Phenotype Ontology",          "both",  None),
    ("genes_to_phenotype.txt",    "HPO gene-to-phenotype annotations",     "human", None),
    ("HP2MP.tsv",                 "HP→MP anchor mapping from MGI",         "human", None),
    ("MGI_PhenoGenoMP.rpt",       "MGD genotype-to-MP annotations",        "mouse", None),
    ("gene_info.gz",              "NCBI gene info (symbol + MGI xrefs)",   "mouse", None),
    ("gene2go.gz",                "NCBI Gene-to-GO mappings",              "both",  "go"),
    (
        "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz",
        "GTEx median gene expression (GCT)",
        "human",
        "expression",
    ),
    (
        "BIOGRID-ALL-5.0.256.tab3.txt",
        "BioGRID protein-protein interactions (tab3)",
        "human",
        "ppi",
    ),
    (
        "basic.obo",
        "Uberon anatomical ontology (basic, no imports)",
        "human",
        "uberon",
    ),
    # uberon/basic.obo is OPTIONAL — the pipeline runs in stub-mode without it.
    # It is therefore NOT in this list; a missing OBO file produces a warning,
    # not a hard exit.  Only add it here if you want to enforce full metadata.
]


def _file_needed(species: str, skip_layers: list[str], req_species: str, layer: str | None) -> bool:
    """Return True if the file is required for this run configuration."""
    if req_species == "human" and species == "mouse":
        return False
    if req_species == "mouse" and species == "human":
        return False
    if layer and layer in skip_layers:
        return False
    return True


def _human_size(nbytes: int) -> str:
    """Format byte count as a human-readable string (e.g. 1.2G, 324K)."""
    for unit in ("B", "K", "M", "G", "T"):
        if abs(nbytes) < 1024 or unit == "T":
            if unit == "B":
                return f"{nbytes}{unit}"
            return f"{nbytes:.1f}{unit}" if nbytes % 1 else f"{nbytes:.0f}{unit}"
        nbytes /= 1024
    return f"{nbytes:.1f}T"


def check_data(args: argparse.Namespace, skip_layers: list[str]) -> None:
    """Validate that all required input files exist before running the pipeline.

    Checks are species-aware and layer-aware: files not needed for the
    current ``--species`` / ``--skip-layers`` configuration are skipped.

    When files are missing the function prints a summary and exits with
    code 1, suggesting ``bash download_data.sh`` to fetch them.
    """
    data_dir = Path(args.data_dir)

    # Build a map of logical names → actual paths for files that may
    # be specified via their own CLI flags (overriding the data_dir default).
    path_overrides: dict[str, Path] = {
        "genes_to_phenotype.txt": Path(args.genes_to_phenotype),
        "MGI_PhenoGenoMP.rpt":    Path(args.mgi),
        "HP2MP.tsv":              Path(args.hp2mp),
        "gene2go.gz":             Path(args.gene2go),
        "gene_info.gz":           Path(args.gene_info),
        "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz": Path(args.gtex),
        "BIOGRID-ALL-5.0.256.tab3.txt":   Path(args.biogrid),
        "basic.obo":       Path(args.uberon_obo),
    }

    print(f"Checking required data files in: {data_dir}/")
    print("─" * 60)

    missing: list[tuple[str, str]] = []
    checked = 0

    for filename, description, req_species, layer in _REQUIRED_FILES:
        if not _file_needed(args.species, skip_layers, req_species, layer):
            continue
        filepath = path_overrides.get(filename, data_dir / filename)
        checked += 1
        if filepath.exists():
            size = _human_size(filepath.stat().st_size)
            print(f"  ✓  {filename:<55s}  {size:>8s}  {description}")
        else:
            print(f"  ✗  {filename:<55s}  {'MISSING':>8s}  {description}")
            missing.append((str(filepath), description))

    print("─" * 60)

    if not missing:
        print(f"All {checked} files present. Ready to run the pipeline.\n")
        return

    print()
    print("=" * 60)
    print("  MISSING DATA FILES")
    print("=" * 60)
    for filepath, description in missing:
        print(f"  ✗  {filepath}")
        print(f"     {description}")
    print()
    print("Fetch all files automatically:")
    print("  bash scripts/download_data.sh")
    print()
    print("Or download individual files manually — see README / run.py docstring.")
    print("=" * 60 + "\n")
    sys.exit(1)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build gene-phenotype BioCypher knowledge graphs."
    )
    # -- Input files -------------------------------------------------------
    p.add_argument("--genes-to-phenotype", default="data/genes_to_phenotype.txt")
    p.add_argument("--mgi",      default="data/MGI_PhenoGenoMP.rpt")
    p.add_argument("--hp2mp",    default="data/HP2MP.tsv",
                   help="HP→MP anchor mapping file from MGI")
    p.add_argument("--gene2go",  default="data/gene2go.gz")
    p.add_argument("--gene-info", default="data/gene_info.gz",
                   help="NCBI gene_info.gz — provides MGI→symbol mapping for mouse "
                        "(https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz)")
    p.add_argument("--data-dir", default="data/",
                   help="Directory containing *.obo files")
    p.add_argument("--gtex", default="data/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz",
                   help="GTEx median TPM GCT file "
                        "(https://gtexportal.org/home/datasets)")
    p.add_argument("--tpm-threshold", type=float, default=5.0,
                   help="Minimum median TPM for a gene-tissue edge (default: 5.0)")
    p.add_argument("--biogrid", default="data/BIOGRID-ALL-5.0.256.tab3.txt",
                   help="BioGRID tab3 file (full dump or human-specific). "
                        "Download from https://downloads.thebiogrid.org/BioGRID")
    p.add_argument("--uberon-obo", default="data/basic.obo",
                   help="Path to uberon/basic.obo "
                        "(download: wget http://purl.obolibrary.org/obo/uberon/basic.obo). "
                        "Default: data/basic.obo")

    # -- Output ------------------------------------------------------------
    p.add_argument("--tsv-out",   default="graph_tsv/")
    p.add_argument("--graph-out", default="biocypher_out/",
                   help="Root output dir. Subdirs human/ and mouse/ are always created "
                        "(e.g. biocypher_out/human/ for --species human).")

    # -- BioCypher config --------------------------------------------------
    p.add_argument("--schema-config",    default=None,
                   help="Explicit schema YAML (single-species builds only)")
    p.add_argument("--biocypher-config", default=None,
                   help="Path to biocypher_config.yaml")

    # -- Species selection -------------------------------------------------
    p.add_argument("--species", choices=["human", "mouse", "both"], default="both",
                   help="Which graph(s) to build (default: both)")

    # -- Behaviour flags ---------------------------------------------------
    p.add_argument("--include-mortality", action="store_true",
                   help="Keep MP:0010768 (mortality/aging) in the HP→MP mapping "
                        "(excluded by default following Barbitoff et al. 2025)")
    p.add_argument("--skip-export",  action="store_true",
                   help="Skip all TSV export steps (TSVs must already exist)")
    p.add_argument("--skip-layers", nargs="*", default=[],
                   help="Layer names to skip, e.g. --skip-layers go expression")
    p.add_argument("--allow-missing-layers", action="store_true",
                   help="Warn instead of failing when a layer's TSV files are missing")
    p.add_argument("--clean-tsv", action="store_true",
                   help="Delete and recreate --tsv-out before exporting")
    p.add_argument("--clean-graph", action="store_true",
                   help="Delete and recreate --graph-out before building")

    return p.parse_args()


def main() -> None:
    args = parse_args()
    tsv_out = Path(args.tsv_out)
    skip_layers = [s.lower() for s in args.skip_layers]

    # Must run before exports so their log calls reach the file handler.
    # Goes to biocypher-log/ (not biocypher_out/) so it survives --clean-graph.
    setup_logging(out_dir=Path("biocypher-log"))

    if args.skip_export:
        print("Skipping all TSV exports (--skip-export set).\n")
    else:
        # Validate input data before doing any work
        check_data(args, skip_layers)

        if args.clean_tsv and tsv_out.exists():
            print(f"Cleaning TSV output directory: {tsv_out.resolve()}")
            shutil.rmtree(tsv_out)
        tsv_out.mkdir(parents=True, exist_ok=True)

        # -- Step 0: HP→MP system-level mapping (human only) ---------------
        # Produces edge_hp_to_mp_top.tsv consumed by Step 1 — must run first.
        if args.species in ("human", "both"):
            print("=" * 60)
            print("Step 0 -- HP→MP system-level mapping")
            print("=" * 60)
            run_mapping(
                genes_to_phenotype_path=args.genes_to_phenotype,
                hp2mp_path=args.hp2mp,
                data_dir=args.data_dir,
                out_dir=tsv_out,
                exclude_mortality=not args.include_mortality,
            )
            print(f"-> {tsv_out.resolve()}\n")

        # -- Step 1: phenotype edges (HPO + MGD) ---------------------------
        print("=" * 60)
        print("Step 1 -- Phenotype edges (HPO + MGD)")
        print("=" * 60)
        run_phenotype_pipeline(
            genes_to_phenotype_path=args.genes_to_phenotype,
            mgi_path=args.mgi,
            data_dir=args.data_dir,
            mapping_dir=tsv_out,
            gene_info_path=args.gene_info,
            out_dir=tsv_out,
            species=args.species,
        )
        print(f"-> {tsv_out.resolve()}\n")

        # -- Step 2: GO annotation edges — human ---------------------------
        if "go" not in skip_layers and args.species in ("human", "both"):
            print("=" * 60)
            print("Step 2 -- GO annotation edges, human (taxon 9606)")
            print("=" * 60)
            run_go_pipeline(
                gene2go_path=args.gene2go,
                genes_to_phenotype_path=args.genes_to_phenotype,
                out_dir=tsv_out,
            )
            print(f"-> {tsv_out.resolve()}\n")

        # -- Step 3: GO annotation edges — mouse ---------------------------
        if "go" not in skip_layers and args.species in ("mouse", "both"):
            print("=" * 60)
            print("Step 3 -- GO annotation edges, mouse (taxon 10090)")
            print("=" * 60)
            run_mouse_go_pipeline(
                gene2go_path=args.gene2go,
                gene_info_path=args.gene_info,
                out_dir=tsv_out,
            )
            print(f"-> {tsv_out.resolve()}\n")

        # -- Step 4: GTEx expression edges — human only --------------------
        if "expression" not in skip_layers and args.species in ("human", "both"):
            print("=" * 60)
            print("Step 4 -- GTEx expression edges, human")
            print("=" * 60)
            run_expression_pipeline(
                gtex_path=args.gtex,
                out_dir=tsv_out,
                tpm_threshold=args.tpm_threshold,
            )
            print(f"-> {tsv_out.resolve()}\n")

        # -- Step 5: BioGRID PPI edges — human only -----------------------
        if "ppi" not in skip_layers and args.species in ("human", "both"):
            print("=" * 60)
            print("Step 5 -- BioGRID protein-protein interactions, human")
            print("=" * 60)
            run_ppi_pipeline(
                biogrid_path=args.biogrid,
                out_dir=tsv_out,
            )
            print(f"-> {tsv_out.resolve()}\n")

        # -- Step 6: Uberon anatomical ontology — human only --------------
        # --uberon-obo is optional; without it, nodes get empty metadata (stub mode).
        # Expression TSV from Step 4 restricts output to tissues in the expression layer.
        if "uberon" not in skip_layers and args.species in ("human", "both"):
            print("=" * 60)
            print("Step 6 -- Uberon anatomical ontology, human")
            print("=" * 60)
            expression_tsv = tsv_out / "edge_human_gene_expressed_in_tissue.tsv"
            run_uberon_pipeline(
                out_dir=tsv_out,
                uberon_obo_path=args.uberon_obo,
                expression_tsv=expression_tsv if expression_tsv.exists() else None,
            )
            print(f"-> {tsv_out.resolve()}\n")

        # -- Step N: add new export steps here -----------------------------

    # -- Final: build BioCypher graph(s) -----------------------------------
    graph_out = Path(args.graph_out)
    if args.clean_graph and graph_out.exists():
        print(f"Cleaning graph output directory: {graph_out.resolve()}")
        shutil.rmtree(graph_out)

    build(
        data_dir=tsv_out,
        out_dir=graph_out,
        species=args.species,
        schema_config_path=args.schema_config,
        biocypher_config_path=args.biocypher_config,
        skip_layers=skip_layers,
        allow_missing_layers=args.allow_missing_layers,
    )


if __name__ == "__main__":
    main()