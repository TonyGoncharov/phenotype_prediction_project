#!/usr/bin/env bash
# download_data.sh — download all required data files for the
# gene-phenotype knowledge graph pipeline.
#
# Usage:
#   bash download_data.sh              # download to default data/ directory
#   bash download_data.sh /path/to/data
#   bash download_data.sh --force      # re-download even if files exist
#
# Requirements: wget
#
# Data sources & licences
# ───────────────────────
#   genes_to_phenotype.txt   HPO (hpo.jax.org)          — see HPO licence
#   MGI_PhenoGenoMP.rpt      MGI/Jackson Lab            — see MGI terms of use
#   mp_hp_mgi_all.sssom.tsv  mapping-commons (GitHub)   — CC BY 4.0
#   gene2go.gz               NCBI FTP                   — public domain
#   gene_info.gz             NCBI FTP                   — public domain
#   hp.obo                   OBO Foundry                — CC BY 4.0
#   mp.obo                   OBO Foundry                — CC BY 4.0
#   go-basic.obo             Gene Ontology Consortium   — CC BY 4.0

set -euo pipefail

# ── Parse arguments ────────────────────────────────────────────────
FORCE=false
DATA_DIR="data"

for arg in "$@"; do
  case "$arg" in
    --force)  FORCE=true ;;
    *)        DATA_DIR="$arg" ;;
  esac
done

mkdir -p "$DATA_DIR"

# ── Colours ────────────────────────────────────────────────────────
if [[ -t 1 ]]; then
  GREEN='\033[0;32m'; YELLOW='\033[0;33m'; CYAN='\033[0;36m'; NC='\033[0m'
else
  GREEN=''; YELLOW=''; CYAN=''; NC=''
fi

# ── Download helper ────────────────────────────────────────────────
# Usage: fetch <URL> <output_filename> [description]
fetch() {
  local url="$1"
  local filename="$2"
  local desc="${3:-$filename}"
  local dest="${DATA_DIR}/${filename}"

  if [[ -f "$dest" && "$FORCE" != true ]]; then
    printf "${YELLOW}  SKIP${NC}  %-30s  (already exists, use --force to re-download)\n" "$filename"
    return 0
  fi

  printf "${CYAN}  GET ${NC}  %-30s  %s\n" "$filename" "$desc"
  wget -q --show-progress -O "$dest" "$url" || {
    echo "  ERROR: failed to download $filename from $url" >&2
    rm -f "$dest"
    return 1
  }
  printf "${GREEN}  OK  ${NC}  %-30s\n" "$filename"
}

echo "Downloading data files to: ${DATA_DIR}/"
echo "════════════════════════════════════════════════"

# ── 1. HPO gene-to-phenotype annotations ──────────────────────────
fetch \
  "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt" \
  "genes_to_phenotype.txt" \
  "HPO gene-to-phenotype annotations"

# ── 2. MGD genotype-to-MP annotations ─────────────────────────────
fetch \
  "https://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt" \
  "MGI_PhenoGenoMP.rpt" \
  "MGD genotype-to-MP annotations"

# ── 3. HP↔MP SSSOM cross-species mapping ──────────────────────────
# NOTE: verify this URL points to the correct file in the
# mapping-commons/mh_mapping_initiative repo. The filename and branch
# may change across releases — check the repo if download fails:
# https://github.com/mapping-commons/mh_mapping_initiative
fetch \
  "https://raw.githubusercontent.com/mapping-commons/mh_mapping_initiative/master/mappings/mp_hp_mgi_all.sssom.tsv" \
  "mp_hp_mgi_all.sssom.tsv" \
  "HP-to-MP SSSOM cross-species mapping"

# ── 4. NCBI gene2go ───────────────────────────────────────────────
fetch \
  "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz" \
  "gene2go.gz" \
  "NCBI Gene-to-GO mappings (all species)"

# ── 5. NCBI gene_info ─────────────────────────────────────────────
fetch \
  "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz" \
  "gene_info.gz" \
  "NCBI gene info (symbol + MGI cross-refs)"

# ── 6. HP ontology ────────────────────────────────────────────────
fetch \
  "https://purl.obolibrary.org/obo/hp.obo" \
  "hp.obo" \
  "Human Phenotype Ontology"

# ── 7. MP ontology ────────────────────────────────────────────────
fetch \
  "https://purl.obolibrary.org/obo/mp.obo" \
  "mp.obo" \
  "Mammalian Phenotype Ontology"

# ── 8. GO ontology (basic) ────────────────────────────────────────
fetch \
  "https://purl.obolibrary.org/obo/go/go-basic.obo" \
  "go-basic.obo" \
  "Gene Ontology (basic)"

# ── Summary ────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════"
echo "Done. Run  bash check_data.sh  to verify all files."