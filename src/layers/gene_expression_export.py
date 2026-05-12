"""src/layers/gene_expression_export.py — GTEx expression pipeline for the KG.

Input : GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz (GCT format)
Output: edge_human_gene_expressed_in_tissue.tsv, node_tissue_gtex.tsv,
        qc_expression_unexpressed_genes.tsv
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

DEFAULT_TPM_THRESHOLD = 5.0

# GCT fixed column names
_GCT_ID_COL  = "Name"
_GCT_SYM_COL = "Description"


# ── Helpers ───────────────────────────────────────────────────────── #

def _tissue_id(tissue_name: str) -> str:
    """Stable node ID for a GTEx tissue name.

    Examples
    --------
    "Adipose - Subcutaneous"           → "GTEX:Adipose_Subcutaneous"
    "Brain - Amygdala"                 → "GTEX:Brain_Amygdala"
    "Skin - Sun Exposed (Lower leg)"   → "GTEX:Skin_Sun_Exposed_Lower_leg"
    """
    normalized = re.sub(r"[^a-zA-Z0-9]+", "_", tissue_name).strip("_")
    return f"GTEX:{normalized}"


# ── Data reader ───────────────────────────────────────────────────── #

def read_gtex_gct(path: str | Path) -> pd.DataFrame:
    """Parse a GTEx GCT file into a long-format DataFrame.

    Returns columns: gene_symbol, tissue_name, median_tpm
    Strips the Ensembl version suffix from gene IDs (not used downstream).
    Rows with missing or non-numeric TPM values are silently dropped.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"GTEx GCT file not found at {path}.\n"
            "Download from https://gtexportal.org/home/datasets\n"
            "  (GTEx_Analysis_v*_gene_median_tpm.gct.gz)"
        )

    logger.debug("Parsing GCT file: %s", path)
    # skiprows=2 skips the "#1.2" line and the dimensions line
    df = pd.read_csv(path, sep="\t", skiprows=2, dtype=str)

    if _GCT_ID_COL not in df.columns or _GCT_SYM_COL not in df.columns:
        raise ValueError(
            f"Unexpected GCT format — expected columns '{_GCT_ID_COL}' and "
            f"'{_GCT_SYM_COL}', got: {list(df.columns[:4])} ..."
        )

    tissue_cols = [c for c in df.columns if c not in (_GCT_ID_COL, _GCT_SYM_COL)]
    if not tissue_cols:
        raise ValueError("No tissue columns found in GCT file.")
    logger.debug("GCT: %d gene rows, %d tissue columns", len(df), len(tissue_cols))

    df_long = df.melt(
        id_vars=[_GCT_SYM_COL],
        value_vars=tissue_cols,
        var_name="tissue_name",
        value_name="median_tpm",
    )
    df_long = df_long.rename(columns={_GCT_SYM_COL: "gene_symbol"})

    before = len(df_long)
    df_long["median_tpm"] = pd.to_numeric(df_long["median_tpm"], errors="coerce")
    df_long = df_long.dropna(subset=["gene_symbol", "tissue_name", "median_tpm"])
    dropped = before - len(df_long)
    if dropped:
        logger.warning("Dropped %d rows with missing/non-numeric TPM values", dropped)

    return df_long[["gene_symbol", "tissue_name", "median_tpm"]]


# ── End-to-end pipeline ───────────────────────────────────────────── #

def run_expression_pipeline(
    gtex_path: str | Path,
    out_dir: str | Path = "./out",
    tpm_threshold: float = DEFAULT_TPM_THRESHOLD,
) -> dict[str, pd.DataFrame]:
    """Run the GTEx expression pipeline for HUMAN genes.

    Filters gene-tissue pairs to those with median TPM >= *tpm_threshold*,
    writes edge and node TSVs, and returns both DataFrames.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting GTEx expression pipeline")
    logger.info("GTEx file     : %s", gtex_path)
    logger.info("Output dir    : %s", out_dir)
    logger.info("TPM threshold : %.1f", tpm_threshold)

    expr = read_gtex_gct(gtex_path)

    all_symbols = set(expr["gene_symbol"].dropna().unique())
    logger.info("Genes in GTEx   : %d", len(all_symbols))
    logger.info("Tissues in GTEx : %d", expr["tissue_name"].nunique())

    # ── Filter by threshold ───────────────────────────────────────── #
    expressed = expr[expr["median_tpm"] >= tpm_threshold].copy()
    expressed["tissue_id"] = expressed["tissue_name"].map(_tissue_id)

    edge_df = (
        expressed
        .groupby(["gene_symbol", "tissue_id", "tissue_name"], as_index=False)
        ["median_tpm"].max()
    )
    edge_df.to_csv(out_dir / "edge_human_gene_expressed_in_tissue.tsv", sep="\t", index=False)
    logger.info("Written: edge_human_gene_expressed_in_tissue.tsv (%d rows)", len(edge_df))
    _dupes = edge_df.duplicated(subset=["gene_symbol", "tissue_id"]).sum()
    if _dupes:
        logger.warning("Duplicate gene→tissue rows in output: %d", _dupes)

    # ── Tissue node table ─────────────────────────────────────────── #
    tissue_df = (
        expressed[["tissue_id", "tissue_name"]]
        .drop_duplicates()
        .sort_values("tissue_name")
        .reset_index(drop=True)
    )
    tissue_df.to_csv(out_dir / "node_tissue_gtex.tsv", sep="\t", index=False)
    logger.info("Written: node_tissue_gtex.tsv (%d rows)", len(tissue_df))

    # ── QC: genes with no tissue above threshold ──────────────────── #
    expressed_symbols = set(edge_df["gene_symbol"].dropna())
    unexpressed = sorted(all_symbols - expressed_symbols)
    pd.DataFrame({"gene_symbol": unexpressed}).to_csv(
        out_dir / "qc_expression_unexpressed_genes.tsv", sep="\t", index=False
    )

    logger.info("Gene→tissue edges      : %d", len(edge_df))
    logger.info("Unique genes expressed : %d", edge_df["gene_symbol"].nunique())
    logger.info("Unique tissues         : %d", edge_df["tissue_id"].nunique())
    logger.info("Genes below threshold  : %d", len(unexpressed))
    if unexpressed:
        logger.debug("Unexpressed genes written to: qc_expression_unexpressed_genes.tsv")

    return {"edges": edge_df, "tissues": tissue_df}


if __name__ == "__main__":
    import argparse
    from src.utils.logging_config import setup_logging

    p = argparse.ArgumentParser(description="Export GTEx gene-expression edges.")
    p.add_argument("--gtex",          required=True, help="GTEx GCT file (.gct or .gct.gz)")
    p.add_argument("--out",           default="./out")
    p.add_argument("--tpm-threshold", type=float, default=DEFAULT_TPM_THRESHOLD,
                   help=f"Minimum median TPM to include an edge (default {DEFAULT_TPM_THRESHOLD})")
    args = p.parse_args()

    setup_logging(out_dir=args.out)
    run_expression_pipeline(
        gtex_path=args.gtex,
        out_dir=args.out,
        tpm_threshold=args.tpm_threshold,
    )
