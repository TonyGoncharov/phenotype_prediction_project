"""export_go_edges.py - GO annotation pipeline for the gene-phenotype KG."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src.utils.gene_info import build_mouse_symbol_maps

logger = logging.getLogger(__name__)

IEA_CODES: frozenset[str] = frozenset({"IEA"})

HUMAN_TAX_ID = "9606"
MOUSE_TAX_ID = "10090"


# ------------------------------------------------------------------ #
# Data readers
# ------------------------------------------------------------------ #

def read_gene2go(path: str | Path, tax_id: str) -> pd.DataFrame:
    """Read NCBI gene2go file, keeping only genes of *tax_id* with non-IEA evidence.

    Reads in 500k-row chunks to avoid loading all ~115M rows into memory at once.
    """
    path = Path(path)
    logger.debug("Reading gene2go from %s (tax_id=%s)", path, tax_id)
    col_names = ["tax_id", "GeneID", "GO_ID", "Evidence", "Qualifier",
                 "GO_term", "PubMed", "Category"]
    chunks: list[pd.DataFrame] = []
    total = 0
    retained = 0
    for chunk in pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        comment="#",
        header=None,
        names=col_names,
        chunksize=500_000,
    ):
        total += len(chunk)
        chunk = chunk[chunk["tax_id"] == tax_id]
        chunk = chunk[~chunk["Evidence"].isin(IEA_CODES)]
        retained += len(chunk)
        if len(chunk):
            chunks.append(chunk)
    df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=col_names)
    logger.debug(
        "gene2go: %d total rows → %d retained after filtering (tax_id=%s, excl. IEA)",
        total, retained, tax_id,
    )

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
    logger.debug("Building NCBI→HGNC map from %s", genes_to_phenotype_path)
    df = pd.read_csv(genes_to_phenotype_path, sep="\t", dtype=str)
    required = {"ncbi_gene_id", "gene_symbol"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"genes_to_phenotype file missing columns for HGNC mapping: {missing}. "
            f"Got: {list(df.columns)}"
        )
    mapping = (
        df[["ncbi_gene_id", "gene_symbol"]]
        .dropna()
        .drop_duplicates()
        .set_index("ncbi_gene_id")["gene_symbol"]
        .to_dict()
    )
    logger.debug("NCBI→HGNC map: %d entries", len(mapping))
    return mapping


# ------------------------------------------------------------------ #
# Edge builders
# ------------------------------------------------------------------ #

def build_gene_go_edges(
    gene2go: pd.DataFrame,
    ncbi_to_symbol: dict[str, str],
    gene_col: str = "gene_symbol",
) -> pd.DataFrame:
    df = gene2go.copy()
    df[gene_col] = df["ncbi_gene_id"].map(ncbi_to_symbol)
    before = len(df)
    df = df.dropna(subset=[gene_col])
    unmapped = before - len(df)
    if unmapped:
        logger.debug("%d gene2go rows dropped: no symbol mapping found", unmapped)
    out_cols = [gene_col, "go_id", "evidence_code", "aspect"]
    if "go_name" in df.columns:
        out_cols.insert(2, "go_name")
    result = df[out_cols].drop_duplicates()
    logger.debug("Gene→GO edges built: %d", len(result))
    return result


# ------------------------------------------------------------------ #
# End-to-end pipelines
# ------------------------------------------------------------------ #

def run_go_pipeline(
    gene2go_path: str | Path,
    genes_to_phenotype_path: str | Path,
    out_dir: str | Path = "./out",
) -> dict[str, pd.DataFrame]:
    """Run the GO annotation pipeline for HUMAN genes (tax_id=9606)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting GO pipeline (human)")
    logger.info("gene2go                : %s", gene2go_path)
    logger.info("genes_to_phenotype     : %s", genes_to_phenotype_path)
    logger.info("Output dir             : %s", out_dir)

    gene2go = read_gene2go(gene2go_path, tax_id=HUMAN_TAX_ID)
    ncbi_to_hgnc = build_ncbi_to_hgnc_map(genes_to_phenotype_path)

    gene_go = build_gene_go_edges(gene2go, ncbi_to_hgnc, gene_col="gene_symbol")
    gene_go.to_csv(out_dir / "edge_human_gene_has_go.tsv", sep="\t", index=False)
    logger.debug("Written: edge_human_gene_has_go.tsv")

    annotated = set(gene_go["gene_symbol"].dropna())
    all_genes = set(ncbi_to_hgnc.values())
    unmapped = sorted(all_genes - annotated)
    pd.DataFrame({"gene_symbol": unmapped}).to_csv(
        out_dir / "qc_go_unmapped_genes.tsv", sep="\t", index=False
    )

    logger.info("Gene→GO edges          : %d", len(gene_go))
    logger.info("Unique genes annotated : %d", gene_go["gene_symbol"].nunique())
    logger.info("Unique GO terms        : %d", gene_go["go_id"].nunique())
    logger.info("Genes without GO annot.: %d", len(unmapped))

    return {"gene_go": gene_go}


def run_mouse_go_pipeline(
    gene2go_path: str | Path,
    gene_info_path: str | Path,
    out_dir: str | Path = "./out",
) -> dict[str, pd.DataFrame]:
    """Run the GO annotation pipeline for MOUSE genes (tax_id=10090)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting GO pipeline (mouse)")
    logger.info("gene2go    : %s", gene2go_path)
    logger.info("gene_info  : %s", gene_info_path)
    logger.info("Output dir : %s", out_dir)

    gene2go = read_gene2go(gene2go_path, tax_id=MOUSE_TAX_ID)

    logger.debug("Building NCBI→mouse symbol map")
    ncbi_to_symbol, _ = build_mouse_symbol_maps(gene_info_path)
    logger.debug("NCBI→mouse symbol map: %d entries", len(ncbi_to_symbol))

    gene_go = build_gene_go_edges(gene2go, ncbi_to_symbol, gene_col="gene_symbol")
    gene_go.to_csv(out_dir / "edge_mouse_gene_has_go.tsv", sep="\t", index=False)
    logger.debug("Written: edge_mouse_gene_has_go.tsv")

    annotated = set(gene_go["gene_symbol"].dropna())
    all_symbols = set(ncbi_to_symbol.values())
    unmapped = sorted(all_symbols - annotated)
    pd.DataFrame({"gene_symbol": unmapped}).to_csv(
        out_dir / "qc_go_unmapped_mouse_genes.tsv", sep="\t", index=False
    )

    logger.info("Mouse gene→GO edges      : %d", len(gene_go))
    logger.info("Unique genes annotated   : %d", gene_go["gene_symbol"].nunique())
    logger.info("Unique GO terms          : %d", gene_go["go_id"].nunique())
    logger.info("Genes without GO annot.  : %d", len(unmapped))

    return {"gene_go": gene_go}


if __name__ == "__main__":
    import argparse
    from src.utils.logging_config import setup_logging

    p = argparse.ArgumentParser(description="Export GO annotation edges.")
    p.add_argument("--gene2go", required=True)
    p.add_argument("--genes-to-phenotype", default=None,
                   help="HPO genes_to_phenotype.txt (required for --species human)")
    p.add_argument("--gene-info", default=None,
                   help="NCBI gene_info.gz (required for --species mouse)")
    p.add_argument("--out", default="./out")
    p.add_argument("--species", choices=["human", "mouse", "both"], default="both")
    args = p.parse_args()

    setup_logging(out_dir=args.out)

    if args.species in ("human", "both"):
        if not args.genes_to_phenotype:
            p.error("--genes-to-phenotype required for human GO export")
        run_go_pipeline(
            gene2go_path=args.gene2go,
            genes_to_phenotype_path=args.genes_to_phenotype,
            out_dir=args.out,
        )

    if args.species in ("mouse", "both"):
        if not args.gene_info:
            p.error("--gene-info required for mouse GO export")
        run_mouse_go_pipeline(
            gene2go_path=args.gene2go,
            gene_info_path=args.gene_info,
            out_dir=args.out,
        )
