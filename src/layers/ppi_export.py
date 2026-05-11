"""
ppi_export.py – BioGRID protein-protein interaction export pipeline (human).

Input : BioGRID tab3 file (BIOGRID-ALL-*.tab3.txt or organism-specific variant)
Output: node_human_protein.tsv, edge_human_gene_encodes_protein.tsv,
        edge_human_ppi.tsv, qc_ppi_unmapped_symbols.tsv
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

HUMAN_TAX_ID = "9606"
PROTEIN_ID_PREFIX = "PROTEIN:"

# BioGRID tab3 column names — two generations of naming are supported.
#
# v3 / v4  (up to ~4.x):  "Official Symbol Interactor A",
#                          "Organism ID Interactor A", "Publication Source", …
# v5       (5.0+):         "Gene Name Interactor A",
#                          "Organism Interactor A", "NCBI Publication ID", …
#
# Each entry maps  logical_field → (v4_name, v5_name).
# read_biogrid() auto-detects which generation the file uses and builds
# the final  {raw_column → internal_name}  rename dict accordingly.
_BIOGRID_FIELD_ALIASES: dict[str, tuple[str, str]] = {
    # logical name            v4 column name                  v5 column name
    "symbol_a":             ("Official Symbol Interactor A",  "Gene Name Interactor A"),
    "symbol_b":             ("Official Symbol Interactor B",  "Gene Name Interactor B"),
    "experimental_system":  ("Experimental System",           "Experimental System"),
    "experimental_system_type": ("Experimental System Type",  "Experimental System Type"),
    "pubmed_id":            ("Publication Source",            "NCBI Publication ID"),
    "tax_id_a":             ("Organism ID Interactor A",      "Organism Interactor A"),
    "tax_id_b":             ("Organism ID Interactor B",      "Organism Interactor B"),
    "throughput":           ("Throughput",                    "Throughput"),
    "source_db":            ("Source Database",               "Source Database"),
}


def _detect_column_map(columns: list[str]) -> dict[str, str]:
    """Return {raw_column → internal_name} for whichever BioGRID generation matches.

    Tries v4 names first, then v5.  Raises ValueError with a diagnostic
    showing the actual columns if neither generation matches.
    """
    col_set = set(columns)

    for gen_idx in (0, 1):          # 0 = v4 names, 1 = v5 names
        rename: dict[str, str] = {}
        for internal, aliases in _BIOGRID_FIELD_ALIASES.items():
            raw = aliases[gen_idx]
            if raw in col_set:
                rename[raw] = internal

        if len(rename) == len(_BIOGRID_FIELD_ALIASES):
            gen_label = "v4" if gen_idx == 0 else "v5"
            logger.debug("BioGRID column format detected: %s", gen_label)
            return rename

    # Neither matched — collect which fields are unresolvable and report
    missing_fields: dict[str, list[str]] = {}
    for internal, aliases in _BIOGRID_FIELD_ALIASES.items():
        if aliases[0] not in col_set and aliases[1] not in col_set:
            missing_fields[internal] = list(aliases)

    raise ValueError(
        f"Cannot map {len(missing_fields)} required BioGRID fields to any known column name.\n\n"
        "Missing fields and their expected names:\n"
        + "\n".join(
            f"  {field}:\n    v4 → '{v4}'\n    v5 → '{v5}'"
            for field, (v4, v5) in _BIOGRID_FIELD_ALIASES.items()
            if field in missing_fields
        )
        + f"\n\nActual columns in file ({len(columns)}):\n  "
        + "\n  ".join(columns[:40])
        + ("\n  …" if len(columns) > 40 else "")
        + "\n\nMake sure you are using a BioGRID tab3-format file."
    )


# ─────────────────────────────────────────────────────────────────────────────
# Data readers
# ─────────────────────────────────────────────────────────────────────────────

def read_biogrid(path: str | Path) -> pd.DataFrame:
    """Read a BioGRID tab3 file (v4 or v5), keeping only the columns we need.

    Accepts both the full ALL-organisms dump and organism-specific files.
    Auto-detects column naming generation (v4 vs v5).
    Returns a DataFrame with normalised internal column names; organism
    filtering is done by the caller so the reader stays reusable.
    """
    path = Path(path)
    logger.info("Reading BioGRID from %s", path)

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        low_memory=False,
    )
    # BioGRID tab3 headers start with '#BioGRID Interaction ID'.
    # Strip the leading '#' so column names are addressable by name.
    df.columns = [c.lstrip("#") for c in df.columns]
    logger.debug("BioGRID raw rows: %d, columns: %d", len(df), len(df.columns))

    rename = _detect_column_map(list(df.columns))

    # In v5 the organism columns contain the full taxonomy string
    # e.g. "9606" or "Homo sapiens (9606)" — normalise to bare NCBI tax ID
    df = df[list(rename)].rename(columns=rename)
    for col in ("tax_id_a", "tax_id_b"):
        if col in df.columns:
            # Extract the numeric part if the cell looks like "Homo sapiens (9606)"
            df[col] = (
                df[col]
                .str.extract(r"\((\d+)\)", expand=False)
                .fillna(df[col])
                .str.strip()
            )

    logger.debug("BioGRID rows loaded: %d", len(df))
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Filtering helpers
# ─────────────────────────────────────────────────────────────────────────────

def filter_human_physical(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only intra-human physical interactions.

    Removes:
    * rows where either interactor is not human (tax_id ≠ 9606)
    * genetic interactions (experimental_system_type == "genetic")
    * self-interactions (symbol_a == symbol_b)
    * rows with missing gene symbols
    """
    before = len(df)
    df = df[
        (df["tax_id_a"] == HUMAN_TAX_ID) &
        (df["tax_id_b"] == HUMAN_TAX_ID)
    ].copy()
    logger.debug("After human-only filter: %d / %d rows", len(df), before)

    df = df[df["experimental_system_type"].str.lower() == "physical"].copy()
    logger.debug("After physical-only filter: %d rows", len(df))

    df = df.dropna(subset=["symbol_a", "symbol_b"])
    df = df[df["symbol_a"] != df["symbol_b"]]
    logger.debug("After drop-NA/self-interactions: %d rows", len(df))
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Node builders
# ─────────────────────────────────────────────────────────────────────────────

def build_protein_nodes(ppi: pd.DataFrame) -> pd.DataFrame:
    symbols = pd.concat([ppi["symbol_a"], ppi["symbol_b"]]).dropna().unique()
    logger.debug("Unique protein symbols: %d", len(symbols))
    return pd.DataFrame({
        "protein_id": [f"{PROTEIN_ID_PREFIX}{s}" for s in symbols],
        "gene_symbol": list(symbols),
        "organism": "Homo sapiens",
    })


def build_gene_encodes_protein(protein_nodes: pd.DataFrame) -> pd.DataFrame:
    """One gene → protein edge per protein node (1-to-1 for human)."""
    return pd.DataFrame({
        "gene_symbol": protein_nodes["gene_symbol"],
        "protein_id":  protein_nodes["protein_id"],
    })


# ─────────────────────────────────────────────────────────────────────────────
# Edge builder
# ─────────────────────────────────────────────────────────────────────────────

def build_ppi_edges(ppi: pd.DataFrame) -> pd.DataFrame:
    """Aggregate all evidence for each unique (protein_a, protein_b) pair.

    The pair is canonicalised (sorted alphabetically) so A–B and B–A are one edge.
    throughput collapses to "Low Throughput" if any low-throughput record exists.
    """
    df = ppi.copy()
    df["protein_id_a"] = PROTEIN_ID_PREFIX + df["symbol_a"]
    df["protein_id_b"] = PROTEIN_ID_PREFIX + df["symbol_b"]

    # Canonicalise pair order so (A,B) == (B,A)
    mask = df["protein_id_a"] > df["protein_id_b"]
    df.loc[mask, ["protein_id_a", "protein_id_b"]] = (
        df.loc[mask, ["protein_id_b", "protein_id_a"]].values
    )

    def _sorted_unique(series: pd.Series) -> list[str]:
        return sorted(series.dropna().unique().tolist())

    def _throughput(series: pd.Series) -> str:
        vals = series.dropna().str.lower()
        return "Low Throughput" if (vals == "low throughput").any() else "High Throughput"

    agg = (
        df.groupby(["protein_id_a", "protein_id_b"])
        .agg(
            experimental_systems      =("experimental_system",      _sorted_unique),
            experimental_system_types =("experimental_system_type", _sorted_unique),
            throughput                =("throughput",                _throughput),
            pubmed_ids                =("pubmed_id",                 _sorted_unique),
            source_db                 =("source_db",                 "first"),
        )
        .reset_index()
    )
    logger.debug("Unique PPI edges after aggregation: %d", len(agg))
    return agg


# ─────────────────────────────────────────────────────────────────────────────
# End-to-end pipeline
# ─────────────────────────────────────────────────────────────────────────────

def run_ppi_pipeline(
    biogrid_path: str | Path,
    out_dir: str | Path = "./out",
) -> dict[str, pd.DataFrame]:
    """Full export pipeline for human PPI data from BioGRID."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("PPI export pipeline (human)")
    logger.info("BioGRID path : %s", biogrid_path)
    logger.info("Output dir   : %s", out_dir)
    logger.info("=" * 60)

    raw   = read_biogrid(biogrid_path)
    human = filter_human_physical(raw)

    if human.empty:
        logger.warning("No human physical interactions found — output files will be empty.")

    protein_nodes        = build_protein_nodes(human)
    gene_encodes_protein = build_gene_encodes_protein(protein_nodes)
    ppi_edges            = build_ppi_edges(human)

    protein_nodes.to_csv(
        out_dir / "node_human_protein.tsv", sep="\t", index=False
    )
    logger.debug("Written: node_human_protein.tsv (%d rows)", len(protein_nodes))

    gene_encodes_protein.to_csv(
        out_dir / "edge_human_gene_encodes_protein.tsv", sep="\t", index=False
    )
    logger.debug("Written: edge_human_gene_encodes_protein.tsv (%d rows)", len(gene_encodes_protein))

    # pipe-separate list columns for TSV compatibility
    ppi_out = ppi_edges.copy()
    for col in ("experimental_systems", "experimental_system_types", "pubmed_ids"):
        ppi_out[col] = ppi_out[col].apply(
            lambda v: "|".join(v) if isinstance(v, list) else v
        )
    ppi_out.to_csv(out_dir / "edge_human_ppi.tsv", sep="\t", index=False)
    logger.debug("Written: edge_human_ppi.tsv (%d rows)", len(ppi_out))

    # QC — shouldn't happen but surfaces deprecated BioGRID symbols
    known_symbols = set(protein_nodes["gene_symbol"])
    all_biogrid_symbols = (
        set(human["symbol_a"].dropna()) | set(human["symbol_b"].dropna())
    )
    unmapped = sorted(all_biogrid_symbols - known_symbols)
    pd.DataFrame({"gene_symbol": unmapped}).to_csv(
        out_dir / "qc_ppi_unmapped_symbols.tsv", sep="\t", index=False
    )

    logger.info("Protein nodes              : %d", len(protein_nodes))
    logger.info("Gene→Protein edges         : %d", len(gene_encodes_protein))
    logger.info("Protein↔Protein edges      : %d", len(ppi_edges))
    logger.info("Unique experimental systems: %s",
                sorted({s for row in ppi_edges["experimental_systems"] for s in row}))
    logger.info("Unmapped BioGRID symbols   : %d", len(unmapped))

    return {
        "protein_nodes":        protein_nodes,
        "gene_encodes_protein": gene_encodes_protein,
        "ppi_edges":            ppi_edges,
    }


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    from src.utils.logging_config import setup_logging

    p = argparse.ArgumentParser(
        description="Export human PPI edges from BioGRID tab3 file."
    )
    p.add_argument(
        "--biogrid",
        required=True,
        metavar="FILE",
        help="Path to BIOGRID-*.tab3.txt (full dump or human-specific).",
    )
    p.add_argument("--out", default="./out", metavar="DIR")
    args = p.parse_args()

    setup_logging(out_dir=args.out)
    run_ppi_pipeline(biogrid_path=args.biogrid, out_dir=args.out)