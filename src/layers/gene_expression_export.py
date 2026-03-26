"""src/layers/gene_expression_export.py — GTEx expression pipeline for the KG.

Data sources
------------
  gtex_median_tpm.gct.gz  - GTEx median gene expression per tissue (GCT format)
                            https://gtexportal.org/home/datasets
                            File: GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz
                            (or whichever GTEx release is current)

GCT format
----------
  Line 1 : #1.2
  Line 2 : <n_genes>\\t<n_tissues>
  Line 3+: Name\\tDescription\\t<tissue_1>\\t<tissue_2>\\t...
            ENSG...\\t<symbol>\\t<tpm>\\t<tpm>\\t...

  Name        — Ensembl gene ID with version suffix (e.g. ENSG00000223972.5)
  Description — HGNC gene symbol

Outputs (written to out_dir)
----------------------------
  edge_human_gene_expressed_in_tissue.tsv
      gene_symbol, tissue_id, tissue_name, median_tpm
      One row per (gene, tissue) pair where median TPM >= threshold.

  node_tissue_gtex.tsv
      tissue_id, tissue_name
      One row per unique tissue (used by adapter for node properties).

  qc_expression_unexpressed_genes.tsv
      gene_symbol
      Genes present in GTEx but with no tissue above the TPM threshold.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

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

    # skiprows=2 skips the "#1.2" line and the dimensions line
    df = pd.read_csv(path, sep="\t", skiprows=2, dtype=str)

    # The first two columns are always Name (Ensembl) and Description (symbol)
    if _GCT_ID_COL not in df.columns or _GCT_SYM_COL not in df.columns:
        raise ValueError(
            f"Unexpected GCT format — expected columns '{_GCT_ID_COL}' and "
            f"'{_GCT_SYM_COL}', got: {list(df.columns[:4])} ..."
        )

    tissue_cols = [c for c in df.columns if c not in (_GCT_ID_COL, _GCT_SYM_COL)]
    if not tissue_cols:
        raise ValueError("No tissue columns found in GCT file.")

    # Wide → long
    df_long = df.melt(
        id_vars=[_GCT_SYM_COL],
        value_vars=tissue_cols,
        var_name="tissue_name",
        value_name="median_tpm",
    )
    df_long = df_long.rename(columns={_GCT_SYM_COL: "gene_symbol"})

    df_long["median_tpm"] = pd.to_numeric(df_long["median_tpm"], errors="coerce")
    df_long = df_long.dropna(subset=["gene_symbol", "tissue_name", "median_tpm"])

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

    Parameters
    ----------
    gtex_path:
        Path to the GTEx GCT file (plain or gzip).
    out_dir:
        Directory to write output TSVs.
    tpm_threshold:
        Minimum median TPM to include a gene-tissue edge (default 5.0).

    Outputs
    -------
    edge_human_gene_expressed_in_tissue.tsv
    node_tissue_gtex.tsv
    qc_expression_unexpressed_genes.tsv
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading GTEx GCT file: {gtex_path}")
    expr = read_gtex_gct(gtex_path)

    all_symbols = set(expr["gene_symbol"].dropna().unique())
    print(f"Genes in GTEx          : {len(all_symbols)}")
    print(f"Tissues in GTEx        : {expr['tissue_name'].nunique()}")

    # ── Filter by threshold ───────────────────────────────────────── #
    expressed = expr[expr["median_tpm"] >= tpm_threshold].copy()
    expressed["tissue_id"] = expressed["tissue_name"].map(_tissue_id)

    # Reorder columns for readability
    edge_df = expressed[["gene_symbol", "tissue_id", "tissue_name", "median_tpm"]]
    edge_df.to_csv(out_dir / "edge_human_gene_expressed_in_tissue.tsv", sep="\t", index=False)

    # ── Tissue node table ─────────────────────────────────────────── #
    tissue_df = (
        expressed[["tissue_id", "tissue_name"]]
        .drop_duplicates()
        .sort_values("tissue_name")
        .reset_index(drop=True)
    )
    tissue_df.to_csv(out_dir / "node_tissue_gtex.tsv", sep="\t", index=False)

    # ── QC: genes with no tissue above threshold ──────────────────── #
    expressed_symbols = set(edge_df["gene_symbol"].dropna())
    unexpressed = sorted(all_symbols - expressed_symbols)
    pd.DataFrame({"gene_symbol": unexpressed}).to_csv(
        out_dir / "qc_expression_unexpressed_genes.tsv", sep="\t", index=False
    )

    print(f"Gene→tissue edges      : {len(edge_df)}")
    print(f"Unique genes expressed : {edge_df['gene_symbol'].nunique()}")
    print(f"Unique tissues         : {edge_df['tissue_id'].nunique()}")
    print(f"TPM threshold          : {tpm_threshold}")
    print(f"Genes below threshold  : {len(unexpressed)}")

    return {"edges": edge_df, "tissues": tissue_df}


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Export GTEx gene-expression edges.")
    p.add_argument("--gtex",          required=True, help="GTEx GCT file (.gct or .gct.gz)")
    p.add_argument("--out",           default="./out")
    p.add_argument("--tpm-threshold", type=float, default=DEFAULT_TPM_THRESHOLD,
                   help=f"Minimum median TPM to include an edge (default {DEFAULT_TPM_THRESHOLD})")
    args = p.parse_args()

    run_expression_pipeline(
        gtex_path=args.gtex,
        out_dir=args.out,
        tpm_threshold=args.tpm_threshold,
    )
