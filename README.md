# Phenotype Prediction Project

A pipeline for building cross-species gene-phenotype knowledge graphs using [BioCypher](https://biocypher.org/). It integrates human and mouse genetic data with phenotype ontologies and Gene Ontology annotations, producing Neo4j-ready import files.

The pipeline builds two parallel graphs - **human** and **mouse** - that share the same Mammalian Phenotype (MP) top-level term nodes, enabling direct cross-species comparison in Neo4j.

## Data sources

| File | Source | Description |
|------|--------|-------------|
| `genes_to_phenotype.txt` | [HPO](https://hpo.jax.org/) | Human gene → HPO phenotype annotations |
| `MGI_PhenoGenoMP.rpt` | [MGI](https://www.informatics.jax.org/) | Mouse genotype → MP annotations |
| `mp_hp_mgi_all.sssom.tsv` | [mapping-commons](https://github.com/mapping-commons/mh_mapping_initiative) | HP ↔ MP cross-species mapping (SSSOM) |
| `gene2go.gz` | [NCBI](https://ftp.ncbi.nlm.nih.gov/gene/DATA/) | Gene → GO term mappings (all species) |
| `gene_info.gz` | [NCBI](https://ftp.ncbi.nlm.nih.gov/gene/DATA/) | Gene symbols and MGI cross-references |
| `hp.obo` | [OBO Foundry](https://purl.obolibrary.org/obo/hp.obo) | Human Phenotype Ontology |
| `mp.obo` | [OBO Foundry](https://purl.obolibrary.org/obo/mp.obo) | Mammalian Phenotype Ontology |
| `go-basic.obo` | [Gene Ontology](https://purl.obolibrary.org/obo/go/go-basic.obo) | Gene Ontology (basic) |

## Graph structure

Each species graph contains three node types and two edge types:

```
[HumanGene / MouseGene]
    ├──(has mp top term)──▶ [MpTopTerm]     (shared across species)
    └──(has go term)──────▶ [GoTerm]
```

Node properties include human-readable names for MP top terms and GO terms, gene symbols, and GO namespace/evidence codes on edges.

## Installation

### Prerequisites

- Python 3.10+
- [uv](https://docs.astral.sh/uv/) (recommended) or pip
- Neo4j 5.x (for graph import)

### Setup

```bash
git clone <repository-url>
cd phenotype_prediction_project
uv sync
```

### Download data

A helper script fetches all required files into `data/`:

```bash
bash src/utils/download_data.sh
```

To verify everything is in place:

```bash
bash src/utils/check_data.sh
```

## Usage

### Full pipeline (both species)

```bash
uv run python run.py
```

This runs these steps:

1. **Phenotype edges** - maps human genes to MP top-level terms via HPO→MP (SSSOM), maps mouse genes to MP top-level terms via MGI.
2. **GO annotations (human)** - links human genes to GO terms from NCBI gene2go (taxon 9606, non-IEA evidence).
3. **GO annotations (mouse)** - links mouse genes to GO terms (taxon 10090, non-IEA evidence).
4. **Graph build** - feeds intermediate TSVs into BioCypher, producing Neo4j-admin-import-ready CSVs.

Planned layers (not yet available):

- **Gene expression data**
- **Protein-protein interaction data**

Output:

```
biocypher_out/
  human/   ← CSVs + neo4j-admin-import-call.sh
  mouse/   ← CSVs + neo4j-admin-import-call.sh
```

### Single species

```bash
uv run python run.py --species human
uv run python run.py --species mouse
```

### Skip layers or steps

```bash
# Build without GO annotations
uv run python run.py --skip-layers go

# Skip TSV export (reuse existing TSVs, rebuild graph only)
uv run python run.py --skip-export --clean-graph
```

### Clean rebuild

```bash
uv run python run.py --clean-tsv --clean-graph
```

### Import into Neo4j

Make sure Neo4j is stopped, then run the generated import script:

```bash
bash biocypher_out/human/neo4j-admin-import-call.sh
```

> **Note:** The path to `neo4j-admin` in the generated script is controlled by `import_call_bin_prefix` in `config/biocypher_config.yaml`. Set it to match your system (e.g. `/opt/homebrew/bin/` on macOS with Homebrew).

### All CLI options

```
uv run python run.py --help
```

| Flag | Description |
|------|-------------|
| `--species {human,mouse,both}` | Which graph(s) to build (default: `both`) |
| `--skip-export` | Skip TSV generation, use existing files |
| `--skip-layers go` | Skip specific layers |
| `--clean-tsv` | Delete and recreate intermediate TSVs |
| `--clean-graph` | Delete and recreate BioCypher output |
| `--all-mappings` | Keep all HP→MP mappings instead of best-ranked only |
| `--data-dir PATH` | Directory containing .obo files (default: `data/`) |
| `--graph-out PATH` | Output directory (default: `biocypher_out/`) |
| `--biocypher-config PATH` | Path to `biocypher_config.yaml` |

## Project structure
```
├── run.py                          # Main entry point
├── config/
│   ├── biocypher_config.yaml       # BioCypher + Neo4j settings
│   ├── schema_config_human.yaml    # Human graph schema (Biolink)
│   └── schema_config_mouse.yaml    # Mouse graph schema (Biolink)
├── src/
│   ├── build_graph.py              # BioCypher graph assembly
│   ├── adapters/
│   │   ├── base.py                 # BaseAdapter with shared I/O
│   │   ├── gene_to_phenotype_adapter.py
│   │   └── gene_ontology_adapter.py
│   ├── layers/
│   │   ├── gene_to_phenotype_export.py  # HPO/MGI → TSV export
│   │   └── gene_ontology_export.py      # gene2go → TSV export
│   └── utils/
│       ├── gene_info.py            # NCBI gene_info parser
│       ├── check_data.sh           # Validate data files
│       └── download_data.sh        # Fetch all data files
├── pyproject.toml                  # Project metadata and dependencies
├── uv.lock                         # Pinned dependency versions
├── data/                           # Input data (not in git, fetched by download_data.sh)
├── LICENSE
└── README.md
```

## Adding a new layer

> **This section is a work in progress.**

<!-- TODO: document the full process with a concrete example -->

In brief, adding a new data layer (e.g. expression, protein interactions) involves four files:

1. **Adapter** - create `src/adapters/<name>_adapter.py` with `Human<Name>Adapter` and `Mouse<Name>Adapter` inheriting from `BaseAdapter`. Set `layer_name` as a class attribute.
2. **Export** - create `src/layers/<name>_export.py` with the pipeline that reads raw data and writes intermediate TSVs.
3. **Registration** - add the adapter classes to `SPECIES_LAYERS` in `src/build_graph.py`.
4. **Schema** - add node/edge definitions to `config/schema_config_human.yaml` and `config/schema_config_mouse.yaml`.

See the existing phenotype and GO layers for reference.