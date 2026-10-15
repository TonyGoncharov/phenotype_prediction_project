# Graph release specification

<!-- Formatting: one sentence per line. See .editorconfig and .markdownlint.json. -->

- **Status:** draft for discussion, contract v0.1.0-draft
- **Date:** 2026-10-15
- **Owners:** `kg` (graph construction) and `ml` (models).
  Changes need a PR (pull request) approved by both sides.

This document defines the format in which `g2p_kg` publishes the graph and `g2p_ml` reads it.
Anything not described here is not part of the contract and must not be relied on.

Keywords: **MUST** / **MUST NOT** are hard requirements; the validator reports a violation as an error.
**SHOULD** is a recommendation; a deviation is a warning.
**MAY** is optional.

## 1. Principles

1. `g2p_kg` and `g2p_ml` interact only through a graph release and the `g2p_contract` package.
   They MUST NOT import each other.
2. Labels (gene -> phenotype class) are stored separately from context edges and MUST NOT appear in `edges`.
3. A gene symbol is never an ID.
4. A release is self-describing: `manifest.json` records everything needed to reproduce it (sources, parameters, code version).
5. A release is deterministic: the same inputs and code produce byte-identical files.
6. Data are stored at the biological level where they are measured (gene, protein, …).
   `kg` does not collapse levels; `ml` may derive simpler views.
7. The format is extensible: new layers are added by the procedure in section 4.4.

## 2. Release layout

One release contains both species, built by the same code with the same parameters.

```text
graph-v0.1.0/
├── human/
│   ├── manifest.json
│   ├── nodes.tsv.gz
│   ├── edges.tsv.gz
│   ├── classes.tsv.gz
│   └── labels.tsv.gz
└── mouse/
    └── (same five files)
```

All `.tsv.gz` files:

- UTF-8, tab-separated, one header row, gzip-compressed;
- no quoting: values MUST NOT contain tabs or newlines;
- an empty string means "no value";
- list values are joined with `|`;
- decimal separator is `.`;
- rows are sorted: `nodes` by `id`; `edges` by (`src`, `relation`, `dst`); `labels` by (`gene_id`, `class_id`); `classes` by `class_id`;
- columns appear in the order given below.

## 3. Nodes: `nodes.tsv.gz`

Contains every node of the context graph **and** every labelled gene, even a gene with no edges.

| Column | Type | Required | Description |
|---|---|---|---|
| `id` | str | yes | unique CURIE (Compact URI) |
| `type` | enum | yes | node type (see below) |
| `name` | str | no | gene symbol or term label |
| `xrefs` | list[str] | no | other IDs, e.g. `HGNC:11998\|Ensembl:ENSG00000141510` |

| `type` | `id` pattern (regex, regular expression) | Example |
|---|---|---|
| `gene` | `^NCBIGene:[1-9][0-9]*$` | `NCBIGene:7157` |
| `protein` | see below | `UniProtKB:P04637` |
| `go_term` | `^GO:[0-9]{7}$` | `GO:0005515` |
| `tissue` | `^UBERON:[0-9]{7}$` | `UBERON:0002107` |

The `protein` pattern follows the UniProt accession format:

```text
^UniProtKB:([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$
```

New layers may add node types; see section 4.4.

Rules:

- Genes use NCBI Gene IDs in both species.
  HGNC (HUGO Gene Nomenclature Committee; human), MGI (Mouse Genome Informatics; mouse) and Ensembl IDs go to `xrefs`.
- Proteins use UniProtKB accessions of reviewed (Swiss-Prot) entries.
  One node per entry; isoforms are not separate nodes.
  A gene may encode several proteins.
  The gene -> protein mapping source goes to the manifest.
- Tissues use Uberon (Uber-anatomy ontology) terms in both species.
  If several source tissues map to one Uberon term, `kg` MUST merge them into one node and record the merge rule in the manifest.
  Sources: GTEx (Genotype-Tissue Expression project; human) and GXD (Gene Expression Database at MGI; mouse).
  Source tissue names MAY go to `xrefs`.
  Tissues without a Uberon mapping are dropped and counted in `qc`.
- **[OPEN]** GTEx cell lines (`Cells_*`) have no suitable Uberon term.
- **[OPEN]** Protein isoforms as separate nodes (e.g. `UniProtKB:P04637-2`).

## 4. Context edges: `edges.tsv.gz`

Context edges describe what is known about genes and their products.
The set of layers is not fixed: we plan to add more (for example, post-translational modifications).
This section lists the current layers and MUST be updated whenever a layer is added or changed (section 4.4).

### 4.1 Columns

Core columns, the same for every layer:

| Column | Type | Required | Description |
|---|---|---|---|
| `src` | str | yes | node `id` |
| `relation` | enum | yes | see section 4.2 |
| `dst` | str | yes | node `id` |
| `layer` | enum | yes | determined by `relation` |

Attribute columns follow the core columns.
Each belongs to one layer and is empty for all other layers.

| Column | Type | Layer | Description |
|---|---|---|---|
| `evidence_codes` | list[str] | `go` | GO (Gene Ontology) evidence codes: how each annotation was made |
| `aspect` | enum | `go` | `P` (biological process), `F` (molecular function) or `C` (cellular component) |
| `median_tpm` | float | `expression` | median expression level, TPM |
| `n_experiments` | int | `expression` | number of experiments (GXD) |
| `n_publications` | int | `ppi` | number of distinct PubMed IDs |
| `throughput` | enum | `ppi` | `low` or `high` |

### 4.2 Layers and relations

Relation names are the same for both species (no `Human`/`Mouse` prefix).

| `layer` | `relation` | `src` -> `dst` | Directed | Required attributes |
|---|---|---|---|---|
| `go` | `gene_has_go_term` | `gene` -> `go_term` | yes | `evidence_codes`, `aspect` |
| `expression` | `gene_expressed_in_tissue` | `gene` -> `tissue` | yes | `median_tpm` |
| `gene_product` | `gene_encodes_protein` | `gene` -> `protein` | yes | – |
| `ppi` | `protein_interacts_with_protein` | `protein` – `protein` | no | `n_publications` |

Layer content:

- `go`: GO annotations of genes.
- `expression`: tissues where a gene is expressed.
- `gene_product`: links from genes to the proteins they encode.
- `ppi`: PPI (protein–protein interactions).

Data are stored at the level where they are measured.
Interactions and modifications are properties of proteins, not genes, so protein-level layers (`ppi` now, PTMs later) connect proteins.
`ml` MAY derive simpler views, such as gene–gene PPI; `kg` does not collapse levels.

`gene_product` is a structural layer: it links genes to the protein-level layers.
Removing a protein-level layer in an ablation removes only that layer's edges, not `gene_product`.

### 4.3 Rules

- An undirected edge is stored once, with `src` < `dst` in string order.
  `ml` adds the reverse direction if it needs it.
- Self-loops (`src` = `dst`) and duplicate (`src`, `relation`, `dst`) triples MUST NOT occur.
  Several evidence codes for one gene–GO pair are merged into `evidence_codes`.
- GO edges are direct annotations only, not propagated to ancestor terms.
  Annotations with the `NOT` qualifier are excluded, as are the evidence codes listed in the manifest (default: `IEA`, Inferred from Electronic Annotation; `ND`, No biological Data available).
- PPI sources such as BioGRID (Biological General Repository for Interaction Datasets) report interactions between genes.
  `kg` maps both genes to their proteins and records the mapping rule in the manifest.
- **[OPEN]** GO annotations with IMP (Inferred from Mutant Phenotype) or IGI (Inferred from Genetic Interaction) evidence are based on mutant phenotypes and may duplicate the labels.
- **[OPEN]** One expression threshold rule for both species.

### 4.4 Adding or changing a layer

A PR that adds a layer MUST, in the same PR:

1. Add the layer and its relations to section 4.2: `src` and `dst` types, direction, required attributes.
2. Add any new node types to section 3, with an ID pattern.
3. Add any new attribute columns to section 4.1.
   An existing column is never reused with a different meaning.
4. Add the layer's sources and build parameters to the manifest (section 7).
5. Update the validator (section 8) and the toy release (section 10).
6. Bump `contract_version`: minor for a new layer; major if an existing layer changes meaning, IDs or required columns (section 9).

Both `kg` and `ml` approve the PR.

**[OPEN]** Candidate layers:

- `go_hierarchy`: `go_term_is_a_go_term`, `go_term` -> `go_term`.
- `anatomy`: `tissue_part_of_tissue`, `tissue` -> `tissue`.
- `ptm`: post-translational modifications on proteins.
  Needs a design for modified sites, which have no standard CURIE.
  Options: an edge from `protein` to a PSI-MOD (Proteomics Standards Initiative Protein Modification Ontology) term with `residue` and `position` attributes, or separate site nodes.

## 5. Classes: `classes.tsv.gz`

Classes are terms of MP (Mammalian Phenotype Ontology).

| Column | Type | Required | Description |
|---|---|---|---|
| `class_id` | str | yes | `^MP:[0-9]{7}$` |
| `name` | str | yes | label from `mp.obo` (the MP ontology file) |
| `set` | enum | yes | `core` or `supplementary` |
| `n_genes` | int | yes | number of labelled genes in this species |

The `core` set is identical in both species.
Proposed content: 24 direct children of `MP:0000001`.

| `class_id` | Name |
|---|---|
| `MP:0002006` | neoplasm |
| `MP:0003631` | nervous system phenotype |
| `MP:0005367` | renal/urinary system phenotype |
| `MP:0005369` | muscle phenotype |
| `MP:0005370` | liver/biliary system phenotype |
| `MP:0005371` | limbs/digits/tail phenotype |
| `MP:0005375` | adipose tissue phenotype |
| `MP:0005376` | homeostasis/metabolism phenotype |
| `MP:0005377` | hearing/vestibular/ear phenotype |
| `MP:0005378` | growth/size/body region phenotype |
| `MP:0005379` | endocrine/exocrine gland phenotype |
| `MP:0005380` | embryo phenotype |
| `MP:0005381` | digestive/alimentary phenotype |
| `MP:0005382` | craniofacial phenotype |
| `MP:0005384` | cellular phenotype |
| `MP:0005385` | cardiovascular system phenotype |
| `MP:0005386` | behavior/neurological phenotype |
| `MP:0005387` | immune system phenotype |
| `MP:0005388` | respiratory system phenotype |
| `MP:0005389` | reproductive system phenotype |
| `MP:0005390` | skeleton phenotype |
| `MP:0005391` | vision/eye phenotype |
| `MP:0005397` | hematopoietic system phenotype |
| `MP:0010771` | integument phenotype |

- `supplementary`: `MP:0010768` mortality/aging.
- Excluded:
  - `MP:0002873` normal phenotype – not an organ system;
  - `MP:0001186` pigmentation and `MP:0005394` taste/olfaction – mouse only.
- A class enters `core` only if **both** species have at least `core_class_min_genes` labelled genes in it.
  **[OPEN]** Threshold value; proposed: 50.

## 6. Labels: `labels.tsv.gz`

Labels come from two phenotype ontologies: HPO (Human Phenotype Ontology, terms `HP:…`) for human and MP for mouse.

| Column | Type | Required | Description |
|---|---|---|---|
| `gene_id` | str | yes | gene `id` from `nodes` |
| `class_id` | str | yes | `class_id` from `classes` |
| `source_ids` | list[str] | yes | source terms: `HP:…` (human), `MP:…` (mouse) |
| `min_distance` | int | no | HP term -> anchor term distance (human only) |

Meaning:

- A label is an observed association between a gene and an organ system.
- A missing label means "unknown", **not** "no association".
- Duplicate (`gene_id`, `class_id`) pairs MUST NOT occur.

Construction (implemented in `kg`; parameters go to the manifest):

- **Human.** HPO `genes_to_phenotype.txt` -> HP terms -> nearest anchor terms in HP2MP (the MGI table that maps HP terms to MP terms) -> MP classes.
  **[OPEN]** Filter by annotation frequency (`frequency` column).
- **Mouse.** MGI annotations -> MP term -> MP classes via the ontology hierarchy.
  Only single-gene genotypes; no transgenes, no conditional alleles.
  **[OPEN]** Which MGI report to use.
- Genes that cannot be mapped to NCBI Gene are dropped and counted in `qc`.

## 7. Manifest: `manifest.json`

```json
{
  "contract_version": "0.1.0",
  "release": "graph-v0.1.0",
  "species": "human",
  "taxon": "NCBITaxon:9606",
  "created_at": "2026-11-01T12:00:00Z",
  "kg_commit": "<git SHA>",
  "sources": [
    {
      "name": "gene2go",
      "file": "gene2go.gz",
      "version": "<version or download date>",
      "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",
      "sha256": "<sha256>"
    }
  ],
  "parameters": {
    "go_excluded_evidence": ["IEA", "ND"],
    "go_exclude_not_qualifier": true,
    "expression_min_tpm": 5.0,
    "gene_to_protein": "<mapping source and rule>",
    "ppi_min_publications": 1,
    "core_class_min_genes": 50,
    "excluded_classes": ["MP:0002873", "MP:0001186", "MP:0005394"]
  },
  "counts": {
    "nodes": {"gene": 0, "protein": 0, "go_term": 0, "tissue": 0},
    "edges": {
      "gene_has_go_term": 0,
      "gene_expressed_in_tissue": 0,
      "gene_encodes_protein": 0,
      "protein_interacts_with_protein": 0
    },
    "labels": {"core": 0, "supplementary": 0}
  },
  "qc": {
    "labeled_genes": 0,
    "unmapped_source_genes": {"<source>": 0},
    "layer_coverage": {"go": 0.0, "expression": 0.0, "gene_product": 0.0, "ppi": 0.0}
  },
  "files": {
    "nodes.tsv.gz": "<sha256>",
    "edges.tsv.gz": "<sha256>",
    "classes.tsv.gz": "<sha256>",
    "labels.tsv.gz": "<sha256>"
  }
}
```

`files` holds the SHA-256 checksum (sha256) of each data file.

`layer_coverage` is the fraction of genes with `core` labels that have at least one edge in the layer.
For a protein-level layer, a gene counts if it encodes (`gene_encodes_protein`) a protein that has an edge in the layer.

## 8. Validator checks

Both sides run the validator from `g2p_contract`: `kg` before publishing, `ml` when loading.
An error (E) blocks the release; a warning (W) does not.

| # | Check | Level |
|---|---|---|
| V1 | All five files exist; sha256 values match `files` in the manifest | E |
| V2 | Columns, their order and types match this spec | E |
| V3 | `id` values are unique and match the pattern for their `type` (section 3) | E |
| V4 | Both ends of every edge exist in `nodes`, with the types the relation requires | E |
| V5 | `relation` is in the vocabulary; `layer` matches it | E |
| V6 | No self-loops, duplicate triples or duplicate labels | E |
| V7 | No label appears in `edges`; `nodes` has no MP nodes | E |
| V8 | Every `gene_id` in `labels` is in `nodes`; every `class_id` is in `classes` | E |
| V9 | `layer_coverage` > 0 for every layer | E |
| V10 | `layer_coverage` ≥ 0.5 for every layer | W |
| V11 | `counts` match the actual numbers | E |
| V12 | The `core` set is the same in `human` and `mouse` of one release | E |
| V13 | The major version of `contract_version` matches the installed `g2p_contract` | E |
| V14 | Rows are sorted (section 2) | W |

V9 and V12 catch two errors of the legacy pipeline: a disconnected mouse expression layer and different class sets in the two species.

## 9. Versioning

- The **contract** follows semver (semantic versioning) in `contract_version`:
  - major – breaking change (a column renamed or removed, a new ID scheme, a changed field meaning);
  - minor – compatible addition (a new optional column, relation or layer);
  - patch – wording fixes, no format change.
- A **graph release** is named `graph-vX.Y.Z` and tied to a git tag with the same name.
  Any change in data, parameters or `kg` code needs a new release.
- `ml` pins the exact release version and the manifest sha256 values in its config.
  Picking the "latest" release automatically is not allowed.

## 10. Storage and delivery

- Working releases are GitHub Releases of this repository, tag `graph-vX.Y.Z`.
  Assets: `graph-vX.Y.Z-human.tar.gz`, `graph-vX.Y.Z-mouse.tar.gz` and `SHA256SUMS`.
  The release notes list the changes since the previous release.
- The final release for publication is also deposited on Zenodo and gets a DOI.
- Releases are not stored in git.
  Git holds only the toy graph `contract/fixtures/toy/`, which has the same layout and is used by the tests of both packages.

## 11. `g2p_contract` API (planned)

Names are provisional; nothing is implemented yet.

- `load_release(path) -> Release` – reads and validates a release; returns the `nodes`, `edges`, `classes`, `labels` tables and the manifest.
- `write_release(path, nodes, edges, classes, labels, manifest)` – writes files in the section 2 format and fills in `counts`, `files` and `layer_coverage`.
- `validate_release(path) -> Report` – runs the section 8 checks.
- Constants: node types, relation vocabulary, layers, ID patterns.

## 12. Open questions

1. Protein isoforms as separate nodes (section 3).
2. GTEx cell lines without a Uberon term (section 3).
3. One `edges` file for all layers, or one file per layer (`edges/<layer>.tsv.gz`), so that new layers do not add columns to a shared table (section 4).
4. Candidate layers: GO and Uberon hierarchies, PTMs and the design of modified sites (section 4.4).
5. GO annotations with IMP/IGI evidence: exclude them, or publish both variants (section 4.3).
6. One expression threshold rule for both species (section 4.3).
7. Value of `core_class_min_genes` (section 5).
8. Filter by HPO annotation frequency (section 6).
9. Source of single-gene MGI genotypes (section 6).
10. Which other edge attributes and node features `ml` needs.
