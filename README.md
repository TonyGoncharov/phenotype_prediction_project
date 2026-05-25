# Disentangling molecular predictors of genotype-to-phenotype relationships

This repository builds heterogeneous human and mouse knowledge graphs for
gene-to-phenotype prediction, then evaluates which molecular evidence layers
help recover Mammalian Phenotype (MP) associations.

The motivating biological question is whether gene neighbourhood information
from Gene Ontology annotations, expression atlases, protein-protein
interactions, and anatomy mappings can improve phenotype prediction beyond the
observed phenotype-label prevalence. The project compares two modelling
tracks:

| Track | Framework | Entry point | Prediction task |
|---|---|---|---|
| Knowledge graph embeddings | PyKEEN, RotatE | `run_kge.py` | Score missing gene-phenotype triples |
| Graph neural network classification | PyTorch, R-GCN, logistic regression | `run_rgcn.py` | Predict MP top-level classes from gene embeddings |

Graph construction uses BioCypher. Dependency management uses `uv`
exclusively.

## Installation

System requirements:

- Linux or macOS (Windows has not been tested).
- Python 3.11 or newer.
- `uv` for dependency management.
- For smoke tests: 4 CPU threads and at least 16 GB RAM.
- For full CPU runs: 16 or more CPU threads and 64 GB RAM are recommended.

Set up the repository:

```bash
git clone <repo>
cd <repo>
uv sync
```

Download raw data into `data/raw/`:

```bash
uv run python scripts/download_data.py
```

Large raw and processed data live under `data/`, which is gitignored.

## Usage

Build BioCypher graphs:

```bash
uv run python run_kg.py 
```

Use `--species` to select the source files and schema. `--species both` builds
human and mouse graphs with the same timestamp, which is useful before
cross-species comparisons.

Run KGE ablations:

```bash
uv run python run_kge.py --species human --config configs/machine_learning/kge_human.yaml

uv run python run_kge.py --species mouse --config configs/machine_learning/kge_mouse.yaml
```

KGE uses RotatE by default. The full config uses seed 42, 128-dimensional
embeddings, 64 negatives per positive, and an 80,000 target-step budget so
conditions receive comparable optimization effort. `--fast` is only a smoke
test preset with reduced epochs and dimensionality.

Run R-GCN ablations:

```bash
uv run python run_rgcn.py --species human --config configs/machine_learning/rgcn_human.yaml

uv run python run_rgcn.py --species mouse --config configs/machine_learning/rgcn_mouse.yaml
```

R-GCN uses phenotype edges only as labels. They are never included in the
encoder message-passing graph. The default config uses seed 42, 128-dimensional
embeddings, 500 encoder epochs, five negatives per positive edge, and a fixed
80/20 gene-level split reused across ablation conditions.

## Input and Output Examples

Raw inputs are expected in `data/raw/`. The downloader fetches:

- `hp.obo` and `mp.obo`
- `genes_to_phenotype.txt`
- `HP2MP.tsv`
- `MGI_PhenoGenoMP.rpt`
- `gene2go.gz` and `gene_info.gz`
- `GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz`
- `gxdrnaseq.tar.gz` and `gxd_uberon_mapping.tsv`
- `BIOGRID-ALL-5.0.256.tab3.txt`
- `basic.obo`

Example human HPO input:

```text
ncbi_gene_id	gene_symbol	hpo_id	hpo_name	frequency	disease_id
10	NAT2	HP:0000007	Autosomal recessive inheritance	-	OMIM:243400
10	NAT2	HP:0001939	Abnormality of metabolism/homeostasis	-	OMIM:243400
```

Example processed human phenotype edge:

```text
gene_id	gene_symbol	mp_top_id	source_hp_ids	min_hp_to_anchor_dist
HGNC:A2ML1	A2ML1	MP:0005377	HP:0000403	5
HGNC:A2ML1	A2ML1	MP:0005387	HP:0000403	4
```

Example processed GO edge:

```text
gene_id	gene_symbol	go_id	go_name	evidence_code	aspect
HGNC:NAT2	NAT2	GO:0004060	arylamine N-acetyltransferase activity	IBA	F
HGNC:NAT2	NAT2	GO:0005515	protein binding	IPI	F
```

Example R-GCN per-phenotype output:

```text
mp_term_id,ap,n_pos_test,n_test,baseline
MP:0003631,0.8764,812,1049,0.7741
MP:0005382,0.7576,691,1049,0.6587
```

Primary generated outputs:

- `data/processed/<species>/`: normalized TSV nodes and edges.
- `data/biocypher/<species>_<timestamp>/`: BioCypher CSVs for Neo4j import.
- `results/ml/<species>/rgcn/<condition>/`: R-GCN checkpoints, logs, and AUPRC
  tables.
- `results/ml/<species>/kge/<model>/<condition>/`: PyKEEN checkpoints,
  metrics, and prediction artifacts.

## Repository Layout

```text
phenotype_prediction_project/
|-- run_kg.py
|-- run_kge.py
|-- run_rgcn.py
|-- configs/
|   |-- knowledge_graph/
|   `-- machine_learning/
|-- src/
|   |-- knowledge_graph/
|   |-- machine_learning/
|   `-- utils/
|-- scripts/
|-- tests/
|-- docs/
`-- data/
```

## Evaluation Contract

- `HumanGeneHasMpTopTerm` and `MouseGeneHasMpTopTerm` edges are labels, not
  message-passing context, in the R-GCN track.
- The `MultiLabelBinarizer` is fitted once before ablations, and the seed-42
  80/20 gene split is reused unchanged across conditions.
- Genes absent from an ablation graph receive zero embeddings instead of being
  silently dropped.
- Parallel R-GCN runs use spawn-based subprocesses.
- KGE evaluation intentionally uses `filter_known=True`.

## Limitations

- Missing gene-phenotype links are not true negatives.
- Human and mouse phenotype annotations come from different experimental and curation regimes.
- HP-to-MP mapping quality affects cross-species comparability.
- MP top-level classes are broad phenotype-system labels, not fine-grained phenotypes.
- High predictive performance does not by itself imply biological mechanism.

## References

- HGNC: HUGO Gene Nomenclature Committee, https://www.genenames.org/
- GTEx Portal, https://gtexportal.org/
- BioGRID, https://thebiogrid.org/
- Gene Ontology, https://geneontology.org/
- Uberon anatomy ontology, https://obophenotype.github.io/uberon/
- MGI and Mouse Genome Database, https://www.informatics.jax.org/
- Mammalian Phenotype Ontology, https://www.informatics.jax.org/vocab/mp_ontology/
- Human Phenotype Ontology, https://hpo.jax.org/
- BioCypher, https://biocypher.org/
- PyKEEN, https://pykeen.readthedocs.io/
- RotatE: Sun et al., 2019.
- R-GCN: Schlichtkrull et al., 2018.
