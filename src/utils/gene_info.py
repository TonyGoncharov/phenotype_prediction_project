"""src/utils/gene_info.py — shared helpers for NCBI gene_info.gz parsing.

This module is intentionally free of any layer-specific imports so that
both the phenotype and gene_ontology layers can use it without creating
a circular or cross-layer dependency.

Public API
----------
build_mouse_symbol_maps(gene_info_path) -> (entrez_to_symbol, mgi_to_symbol)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

MOUSE_TAX_ID = "10090"


def build_mouse_symbol_maps(
    gene_info_path: str | Path,
) -> tuple[dict[str, str], dict[str, str]]:
    """Parse NCBI gene_info.gz and return two lookup dicts for mouse genes.

    Args:
        gene_info_path: Path to NCBI gene_info.gz
                        (https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz)

    Returns:
        entrez_to_symbol : dict mapping Entrez gene ID  → official symbol
        mgi_to_symbol    : dict mapping MGI accession ID → official symbol
                           (e.g. "MGI:98834" → "Trp53")

    Notes:
        - Only mouse genes (tax_id == 10090) are loaded.
        - The dbXrefs column contains pipe-separated cross-references such as
          "MGI:MGI:98834|Ensembl:...".  Every MGI accession is mapped to the
          gene symbol using vectorised pandas operations (no iterrows).
        - Rows where Symbol equals "-" (NCBI placeholder) are skipped.
        - Duplicate "MGI:MGI:" prefix is normalised to "MGI:".
    """
    gene_info_path = Path(gene_info_path)
    if not gene_info_path.exists():
        raise FileNotFoundError(gene_info_path)

    df = pd.read_csv(
        gene_info_path,
        sep="\t",
        dtype=str,
        usecols=["#tax_id", "GeneID", "Symbol", "dbXrefs"],
    )

    # Mouse only, skip NCBI placeholder rows
    mouse = df[(df["#tax_id"] == MOUSE_TAX_ID) & (df["Symbol"] != "-")].copy()

    entrez_to_symbol: dict[str, str] = dict(zip(mouse["GeneID"], mouse["Symbol"]))

    # ── MGI ID parsing — vectorised, no iterrows ───────────────────── #
    #
    # Strategy:
    #   1. Drop rows with no dbXrefs.
    #   2. Split the pipe-separated xref string into a list column.
    #   3. explode() → one xref per row.
    #   4. Keep only xrefs that start with "MGI:".
    #   5. Normalise "MGI:MGI:XXXXXX" → "MGI:XXXXXX".
    #   6. Build the dict from the resulting (xref, Symbol) pairs.
    #
    mgi_df = (
        mouse[["Symbol", "dbXrefs"]]
        .dropna(subset=["dbXrefs"])
        .assign(xref=lambda d: d["dbXrefs"].str.split("|"))
        .explode("xref")
    )
    mgi_df["xref"] = mgi_df["xref"].str.strip()
    mgi_df = mgi_df[mgi_df["xref"].str.startswith("MGI:", na=False)].copy()

    # Normalise double prefix: "MGI:MGI:98834" → "MGI:98834"
    double = mgi_df["xref"].str.startswith("MGI:MGI:")
    mgi_df.loc[double, "xref"] = "MGI:" + mgi_df.loc[double, "xref"].str[len("MGI:MGI:"):]

    mgi_to_symbol: dict[str, str] = (
        mgi_df[["xref", "Symbol"]]
        .drop_duplicates(subset=["xref"])
        .set_index("xref")["Symbol"]
        .to_dict()
    )

    return entrez_to_symbol, mgi_to_symbol