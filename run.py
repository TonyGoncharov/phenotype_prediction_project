"""run.py - single entry point for the gene-phenotype knowledge graph pipeline.

Usage
-----
  uv run python run.py                          # full pipeline, both graphs
  uv run python run.py --species human          # human graph only
  uv run python run.py --species mouse          # mouse graph only
  uv run python run.py --skip-export            # skip TSV generation
  uv run python run.py --clean-tsv --clean-graph
  uv run python run.py --help

Required data files (put in data/)
-------------------------------------
  genes_to_phenotype.txt   - HPO gene→phenotype annotations
  MGI_PhenoGenoMP.rpt      - MGD genotype→MP annotations
  mp_hp_mgi_all.sssom.tsv  - HP→MP SSSOM cross-species mapping
  gene2go.gz               - NCBI Gene→GO mappings (all species)
  gene_info.gz             - NCBI gene info (symbol + MGI cross-refs)
                             https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
  data/hp.obo              - HP ontology
  data/mp.obo              - MP ontology
  data/go-basic.obo        - GO ontology

Output layout (--species both, the default)
---------------------------------------------
  biocypher_out/
    human/   <- human graph CSVs + neo4j-admin-import-call.sh
    mouse/   <- mouse graph CSVs + neo4j-admin-import-call.sh

Adding a new data layer
-----------------------
  1. Create src/layers/<n>/adapter.py with Human<N>Adapter and Mouse<N>Adapter
     inheriting BaseAdapter.  Set layer_name on each subclass.
  2. Create src/layers/<n>/export.py with the export pipeline.
  3. Register both adapter classes in src/pipeline/build_graph.py → SPECIES_LAYERS.
  4. Add schema sections to config/schema_config_human.yaml / _mouse.yaml.
  5. Add the export call below (Step N).
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from src.layers.gene_to_phenotype_export import run_pipeline as run_phenotype_pipeline
from src.layers.gene_ontology_export import run_go_pipeline, run_mouse_go_pipeline
from src.build_graph import build


# ── Data-file validation ──────────────────────────────────────────── #

# (filename, description, needed_for_species, layer)
# species: "human", "mouse", or "both"
# layer:   None = always needed; "go" = only when GO layer is active
_REQUIRED_FILES: list[tuple[str, str, str, str | None]] = [
    ("hp.obo",                   "Human Phenotype Ontology",            "both",  None),
    ("mp.obo",                   "Mammalian Phenotype Ontology",        "both",  None),
    ("genes_to_phenotype.txt",   "HPO gene-to-phenotype annotations",   "human", None),
    ("mp_hp_mgi_all.sssom.tsv",  "HP-to-MP SSSOM cross-species mapping","human", None),
    ("MGI_PhenoGenoMP.rpt",      "MGD genotype-to-MP annotations",      "mouse", None),
    ("gene_info.gz",             "NCBI gene info (symbol + MGI xrefs)", "mouse", None),
    ("gene2go.gz",               "NCBI Gene-to-GO mappings",            "both",  "go"),
    ("go-basic.obo",             "Gene Ontology (basic)",               "both",  "go"),
]


def _file_needed(species: str, skip_layers: list[str], req_species: str, layer: str | None) -> bool:
    """Return True if the file is required for this run configuration."""
    # Species filter
    if req_species == "human" and species == "mouse":
        return False
    if req_species == "mouse" and species == "human":
        return False
    # Layer filter
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
        "genes_to_phenotype.txt":  Path(args.genes_to_phenotype),
        "MGI_PhenoGenoMP.rpt":     Path(args.mgi),
        "mp_hp_mgi_all.sssom.tsv": Path(args.sssom),
        "gene2go.gz":              Path(args.gene2go),
        "gene_info.gz":            Path(args.gene_info),
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
            print(f"  ✓  {filename:<30s}  {size:>8s}  {description}")
        else:
            print(f"  ✗  {filename:<30s}  {'MISSING':>8s}  {description}")
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
    print("  bash download_data.sh")
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
    p.add_argument("--mgi",                default="data/MGI_PhenoGenoMP.rpt")
    p.add_argument("--sssom",              default="data/mp_hp_mgi_all.sssom.tsv")
    p.add_argument("--gene2go",            default="data/gene2go.gz")
    p.add_argument("--gene-info",          default="data/gene_info.gz",
                   help="NCBI gene_info.gz — provides MGI→symbol mapping for mouse "
                        "(https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz)")
    p.add_argument("--data-dir",           default="data/",
                   help="Directory containing *.obo files")

    # -- Output ------------------------------------------------------------
    p.add_argument("--tsv-out",   default="graph_tsv/")
    p.add_argument("--graph-out", default="biocypher_out/",
                   help="Root output dir. Subdirs human/ and mouse/ are created "
                        "automatically when --species both.")

    # -- BioCypher config --------------------------------------------------
    p.add_argument("--schema-config",    default=None,
                   help="Explicit schema YAML (single-species builds only)")
    p.add_argument("--biocypher-config", default=None,
                   help="Path to biocypher_config.yaml")

    # -- Species selection -------------------------------------------------
    p.add_argument("--species", choices=["human", "mouse", "both"], default="both",
                   help="Which graph(s) to build (default: both)")

    # -- Behaviour flags ---------------------------------------------------
    p.add_argument("--all-mappings", action="store_true",
                   help="Keep all HP->MP mappings instead of best-ranked only")
    p.add_argument("--skip-export",  action="store_true",
                   help="Skip all TSV export steps (TSVs must already exist)")
    p.add_argument("--skip-layers", nargs="*", default=[],
                   help="Layer names to skip, e.g. --skip-layers go")
    p.add_argument("--clean-tsv", action="store_true",
                   help="Delete and recreate --tsv-out before exporting")
    p.add_argument("--clean-graph", action="store_true",
                   help="Delete and recreate --graph-out before building")

    return p.parse_args()


def main() -> None:
    args = parse_args()
    tsv_out = Path(args.tsv_out)
    skip_layers = [s.lower() for s in args.skip_layers]

    if args.skip_export:
        print("Skipping all TSV exports (--skip-export set).\n")
    else:
        # Validate input data before doing any work
        check_data(args, skip_layers)

        if args.clean_tsv and tsv_out.exists():
            print(f"Cleaning TSV output directory: {tsv_out.resolve()}")
            shutil.rmtree(tsv_out)
        tsv_out.mkdir(parents=True, exist_ok=True)

        # -- Step 1: phenotype edges (HPO + MGD) ---------------------------
        # HPO is human-only; MGD is mouse-only — pass species so the
        # pipeline skips whichever half is not needed.
        print("=" * 60)
        print("Step 1 -- Phenotype edges (HPO + MGD)")
        print("=" * 60)
        run_phenotype_pipeline(
            genes_to_phenotype_path=args.genes_to_phenotype,
            mgi_path=args.mgi,
            sssom_path=args.sssom,
            data_dir=args.data_dir,
            gene_info_path=args.gene_info,
            out_dir=tsv_out,
            best_only=not args.all_mappings,
            species=args.species,  # ← new: skip human/mouse halves when not needed
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
                data_dir=args.data_dir,
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
                data_dir=args.data_dir,
                out_dir=tsv_out,
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
    )


if __name__ == "__main__":
    main()