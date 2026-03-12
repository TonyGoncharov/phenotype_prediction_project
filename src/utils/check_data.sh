#!/usr/bin/env bash
# check_data.sh — verify that all required data files are present.
#
# Usage:
#   bash check_data.sh              # checks default data/ directory
#   bash check_data.sh /path/to/data

set -euo pipefail

DATA_DIR="${1:-data}"

# ── Required files ─────────────────────────────────────────────────
# Each entry: "relative/path|description|source"
REQUIRED_FILES=(
  "genes_to_phenotype.txt|HPO gene-to-phenotype annotations|https://github.com/obophenotype/human-phenotype-ontology/releases"
  "MGI_PhenoGenoMP.rpt|MGD genotype-to-MP annotations|https://www.informatics.jax.org/downloads/reports/index.html"
  "mp_hp_mgi_all.sssom.tsv|HP-to-MP SSSOM cross-species mapping|https://github.com/mapping-commons/mh_mapping_initiative"
  "gene2go.gz|NCBI Gene-to-GO mappings (all species)|https://ftp.ncbi.nlm.nih.gov/gene/DATA/"
  "gene_info.gz|NCBI gene info (symbol + MGI cross-refs)|https://ftp.ncbi.nlm.nih.gov/gene/DATA/"
  "hp.obo|Human Phenotype Ontology|https://purl.obolibrary.org/obo/hp.obo"
  "mp.obo|Mammalian Phenotype Ontology|https://purl.obolibrary.org/obo/mp.obo"
  "go-basic.obo|Gene Ontology (basic)|https://purl.obolibrary.org/obo/go/go-basic.obo"
)

# ── Colours (disabled if stdout is not a terminal) ─────────────────
if [[ -t 1 ]]; then
  GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[0;33m'; NC='\033[0m'
else
  GREEN=''; RED=''; YELLOW=''; NC=''
fi

# ── Check loop ─────────────────────────────────────────────────────
echo "Checking required data files in: ${DATA_DIR}/"
echo "────────────────────────────────────────────────"

missing=0
total=0

for entry in "${REQUIRED_FILES[@]}"; do
  IFS='|' read -r filename description source <<< "$entry"
  total=$((total + 1))
  filepath="${DATA_DIR}/${filename}"

  if [[ -f "$filepath" ]]; then
    size=$(du -h "$filepath" | cut -f1 | xargs)
    printf "${GREEN}  ✓${NC}  %-30s  %8s  %s\n" "$filename" "$size" "$description"
  else
    printf "${RED}  ✗${NC}  %-30s  %8s  %s\n" "$filename" "MISSING" "$description"
    printf "${YELLOW}     └─ Download: %s${NC}\n" "$source"
    missing=$((missing + 1))
  fi
done

# ── Summary ────────────────────────────────────────────────────────
echo "────────────────────────────────────────────────"
if [[ $missing -eq 0 ]]; then
  printf "${GREEN}All %d files present. Ready to run the pipeline.${NC}\n" "$total"
  exit 0
else
  printf "${RED}%d of %d files missing.${NC}\n" "$missing" "$total"
  echo "Run  bash download_data.sh  to fetch missing files."
  exit 1
fi