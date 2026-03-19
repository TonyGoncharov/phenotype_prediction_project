"""export_go_edges.py - GO annotation pipeline for the gene-phenotype KG.

Data sources
------------
  gene2go       - NCBI Gene → GO term mappings
                  https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
                  Filtered to: human (tax_id=9606) and/or mouse (tax_id=10090),
                  evidence code != IEA
  gene_info     - NCBI gene info (symbol + MGI cross-refs for mouse)
                  https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
  go-basic.obo  - Gene Ontology (obo format) — used only to validate term IDs
                  https://purl.obolibrary.org/obo/go/go-basic.obo

Outputs (written to out_dir)
----------------------------
  edge_human_gene_has_go.tsv      - HumanGene → GOTerm
  edge_mouse_gene_has_go.tsv      - MouseGene → GOTerm
  qc_go_unmapped_genes.tsv        - human genes with no GO annotation (QC)
  qc_go_unmapped_mouse_genes.tsv  - mouse genes with no GO annotation (QC)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.utils.gene_info import build_mouse_symbol_maps

IEA_CODES: frozenset[str] = frozenset({"IEA"})

HUMAN_TAX_ID = "9606"
MOUSE_TAX_ID = "10090"


# ------------------------------------------------------------------ #
# Ontology helpers
# ------------------------------------------------------------------ #


# ------------------------------------------------------------------ #
# Data readers
# ------------------------------------------------------------------ #

def read_gene2go(path: str | Path, tax_id: str) -> pd.DataFrame:
    """Read NCBI gene2go file, keeping only genes of *tax_id* with non-IEA evidence.

    Returns a DataFrame with columns:
      ncbi_gene_id, go_id, evidence_code, aspect
    where aspect in {P, F, C}.
    """
    path = Path(path)
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        comment="#",
        header=None,
        names=["tax_id", "GeneID", "GO_ID", "Evidence", "Qualifier",
               "GO_term", "PubMed", "Category"],
    )
    df = df[df["tax_id"] == tax_id].copy()
    df = df[~df["Evidence"].isin(IEA_CODES)].copy()
    df = df.rename(columns={
        "GeneID": "ncbi_gene_id",
        "GO_ID": "go_id",
        "Evidence": "evidence_code",
        "Category": "aspect_long",
        "GO_term": "go_name",
    })
    aspect_map = {
        "biological_process": "P",
        "molecular_function": "F",
        "cellular_component": "C",
        "Process": "P",
        "Function": "F",
        "Component": "C",
    }
    df["aspect"] = df["aspect_long"].map(aspect_map).fillna(df["aspect_long"])
    return df[["ncbi_gene_id", "go_id", "go_name", "evidence_code", "aspect"]].drop_duplicates()


def build_ncbi_to_hgnc_map(genes_to_phenotype_path: str | Path) -> dict[str, str]:
    """Build {ncbi_gene_id → hgnc_symbol} from the HPO annotation file."""
    df = pd.read_csv(genes_to_phenotype_path, sep="\t", dtype=str)
    required = {"ncbi_gene_id", "gene_symbol"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"genes_to_phenotype file missing columns for HGNC mapping: {missing}. "
            f"Got: {list(df.columns)}"
        )
    return (
        df[["ncbi_gene_id", "gene_symbol"]]
        .dropna()
        .drop_duplicates()
        .set_index("ncbi_gene_id")["gene_symbol"]
        .to_dict()
    )


# ------------------------------------------------------------------ #
# Edge builders (shared between species)
# ------------------------------------------------------------------ #

def build_gene_go_edges(
    gene2go: pd.DataFrame,
    ncbi_to_symbol: dict[str, str],
    gene_col: str = "gene_symbol",
) -> pd.DataFrame:
    """Map gene2go annotations to gene symbols.

    Returns a DataFrame with columns:
      <gene_col>, go_id, go_name, evidence_code, aspect
    """
    df = gene2go.copy()
    df[gene_col] = df["ncbi_gene_id"].map(ncbi_to_symbol)
    df = df.dropna(subset=[gene_col])
    out_cols = [gene_col, "go_id", "evidence_code", "aspect"]
    if "go_name" in df.columns:
        out_cols.insert(2, "go_name")
    return df[out_cols].drop_duplicates()


# ------------------------------------------------------------------ #
# End-to-end pipelines
# ------------------------------------------------------------------ #

def run_go_pipeline(
    gene2go_path: str | Path,
    genes_to_phenotype_path: str | Path,
    data_dir: str | Path,
    out_dir: str | Path = "./out",
) -> dict[str, pd.DataFrame]:
    """Run the GO annotation pipeline for HUMAN genes (tax_id=9606).

    Outputs written to *out_dir*:
        edge_human_gene_has_go.tsv   - gene_symbol, go_id, evidence_code, aspect
        qc_go_unmapped_genes.tsv     - human gene symbols with no GO annotation
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Reading gene2go annotations (human)...")
    gene2go = read_gene2go(gene2go_path, tax_id=HUMAN_TAX_ID)

    print("Building NCBI → HGNC symbol map...")
    ncbi_to_hgnc = build_ncbi_to_hgnc_map(genes_to_phenotype_path)

    gene_go = build_gene_go_edges(gene2go, ncbi_to_hgnc, gene_col="gene_symbol")
    gene_go.to_csv(out_dir / "edge_human_gene_has_go.tsv", sep="\t", index=False)

    annotated = set(gene_go["gene_symbol"].dropna())
    all_genes = set(ncbi_to_hgnc.values())
    unmapped = sorted(all_genes - annotated)
    pd.DataFrame({"gene_symbol": unmapped}).to_csv(
        out_dir / "qc_go_unmapped_genes.tsv", sep="\t", index=False
    )

    print(f"Gene→GO edges          : {len(gene_go)}")
    print(f"Unique genes annotated : {gene_go['gene_symbol'].nunique()}")
    print(f"Unique GO terms        : {gene_go['go_id'].nunique()}")
    print(f"Genes without GO annot.: {len(unmapped)}")

    return {"gene_go": gene_go}


def run_mouse_go_pipeline(
    gene2go_path: str | Path,
    gene_info_path: str | Path,
    data_dir: str | Path,
    out_dir: str | Path = "./out",
) -> dict[str, pd.DataFrame]:
    """Run the GO annotation pipeline for MOUSE genes (tax_id=10090).

    Outputs written to *out_dir*:
        edge_mouse_gene_has_go.tsv        - gene_symbol, go_id, evidence_code, aspect
        qc_go_unmapped_mouse_genes.tsv    - mouse gene symbols with no GO annotation
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Reading gene2go annotations (mouse)...")
    gene2go = read_gene2go(gene2go_path, tax_id=MOUSE_TAX_ID)

    print("Building NCBI → mouse gene symbol map...")
    ncbi_to_symbol, _ = build_mouse_symbol_maps(gene_info_path)

    gene_go = build_gene_go_edges(gene2go, ncbi_to_symbol, gene_col="gene_symbol")
    gene_go.to_csv(out_dir / "edge_mouse_gene_has_go.tsv", sep="\t", index=False)

    annotated = set(gene_go["gene_symbol"].dropna())
    all_symbols = set(ncbi_to_symbol.values())
    unmapped = sorted(all_symbols - annotated)
    pd.DataFrame({"gene_symbol": unmapped}).to_csv(
        out_dir / "qc_go_unmapped_mouse_genes.tsv", sep="\t", index=False
    )

    print(f"Mouse gene→GO edges      : {len(gene_go)}")
    print(f"Unique genes annotated   : {gene_go['gene_symbol'].nunique()}")
    print(f"Unique GO terms          : {gene_go['go_id'].nunique()}")
    print(f"Genes without GO annot.  : {len(unmapped)}")

    return {"gene_go": gene_go}


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Export GO annotation edges.")
    p.add_argument("--gene2go", required=True)
    p.add_argument("--genes-to-phenotype", default=None,
                   help="HPO genes_to_phenotype.txt (required for --species human)")
    p.add_argument("--gene-info", default=None,
                   help="NCBI gene_info.gz (required for --species mouse)")
    p.add_argument("--data-dir", required=True)
    p.add_argument("--out", default="./out")
    p.add_argument("--species", choices=["human", "mouse", "both"], default="both")
    args = p.parse_args()

    if args.species in ("human", "both"):
        if not args.genes_to_phenotype:
            p.error("--genes-to-phenotype required for human GO export")
        run_go_pipeline(
            gene2go_path=args.gene2go,
            genes_to_phenotype_path=args.genes_to_phenotype,
            data_dir=args.data_dir,
            out_dir=args.out,
        )

    if args.species in ("mouse", "both"):
        if not args.gene_info:
            p.error("--gene-info required for mouse GO export")
        run_mouse_go_pipeline(
            gene2go_path=args.gene2go,
            gene_info_path=args.gene_info,
            data_dir=args.data_dir,
            out_dir=args.out,
        )