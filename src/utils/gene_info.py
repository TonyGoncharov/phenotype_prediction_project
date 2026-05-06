"""src/utils/gene_info.py — shared helpers for NCBI gene_info.gz parsing."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

MOUSE_TAX_ID = "10090"


def build_mouse_symbol_maps(
    gene_info_path: str | Path,
) -> tuple[dict[str, str], dict[str, str]]:
    """Parse NCBI gene_info.gz and return two lookup dicts for mouse genes.

    Returns:
        entrez_to_symbol : dict mapping Entrez gene ID  → official symbol
        mgi_to_symbol    : dict mapping MGI accession ID → official symbol
    """
    gene_info_path = Path(gene_info_path)
    if not gene_info_path.exists():
        raise FileNotFoundError(gene_info_path)

    logger.debug("Reading gene_info from %s", gene_info_path)
    df = pd.read_csv(
        gene_info_path,
        sep="\t",
        dtype=str,
        usecols=["#tax_id", "GeneID", "Symbol", "dbXrefs"],
    )

    mouse = df[(df["#tax_id"] == MOUSE_TAX_ID) & (df["Symbol"] != "-")].copy()
    logger.debug("Mouse gene_info rows: %d", len(mouse))

    entrez_to_symbol: dict[str, str] = dict(zip(mouse["GeneID"], mouse["Symbol"]))

    mgi_df = (
        mouse[["Symbol", "dbXrefs"]]
        .dropna(subset=["dbXrefs"])
        .assign(xref=lambda d: d["dbXrefs"].str.split("|"))
        .explode("xref")
    )
    mgi_df["xref"] = mgi_df["xref"].str.strip()
    mgi_df = mgi_df[mgi_df["xref"].str.startswith("MGI:", na=False)].copy()

    double = mgi_df["xref"].str.startswith("MGI:MGI:")
    mgi_df.loc[double, "xref"] = "MGI:" + mgi_df.loc[double, "xref"].str[len("MGI:MGI:"):]

    mgi_to_symbol: dict[str, str] = (
        mgi_df[["xref", "Symbol"]]
        .drop_duplicates(subset=["xref"])
        .set_index("xref")["Symbol"]
        .to_dict()
    )

    logger.debug(
        "Symbol maps built — Entrez→symbol: %d entries, MGI→symbol: %d entries",
        len(entrez_to_symbol), len(mgi_to_symbol),
    )
    return entrez_to_symbol, mgi_to_symbol
