"""src/adapters/gene_expression_adapter.py — GTEx expression layer adapter."""

from __future__ import annotations

from typing import Generator

from src.adapters.base import BaseAdapter, NodeTuple, EdgeTuple


def _sanitize(value: str) -> str:
    """Remove single quotes that break neo4j-admin CSV import (quote char = ')."""
    return value.replace("'", "")


class HumanExpressionAdapter(BaseAdapter):
    """Adapter for the GTEx gene-expression layer (human only).

    Reads edge_human_gene_expressed_in_tissue.tsv and node_tissue_gtex.tsv
    produced by gene_expression_export.run_expression_pipeline().

    Emits
    -----
    Nodes : human gene, tissue
    Edges : human gene expressed in tissue  (property: median_tpm)
    """

    layer_name = "expression"

    _EDGE_FILE   = "edge_human_gene_expressed_in_tissue.tsv"
    _TISSUE_FILE = "node_tissue_gtex.tsv"
    _GENE_LABEL  = "human gene"
    _EDGE_LABEL  = "human gene expressed in tissue"
    _ID_PREFIX   = "HGNC:"
    _ORGANISM    = "Homo sapiens"

    REQUIRED_FILES: list[str] = [_EDGE_FILE]

    def __init__(self, data_dir):
        super().__init__(data_dir)
        missing = [f for f in self.REQUIRED_FILES if not (self.data_dir / f).exists()]
        if missing:
            raise FileNotFoundError(
                f"Expression layer: missing TSV files in {self.data_dir}:\n  "
                + "\n  ".join(missing)
            )

    # ── Nodes ──────────────────────────────────────────────────────── #

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        yield from self._gene_nodes()
        yield from self._tissue_nodes()

    def _gene_nodes(self) -> Generator[NodeTuple, None, None]:
        df = self._read(self._EDGE_FILE)
        for symbol in df["gene_symbol"].dropna().unique():
            yield (f"{self._ID_PREFIX}{symbol}", self._GENE_LABEL, {
                "symbol": symbol,
                "organism": self._ORGANISM,
            })

    def _tissue_nodes(self) -> Generator[NodeTuple, None, None]:
        """Yield one Tissue node per unique GTEx tissue.

        Falls back to deriving names from the edge file if the separate
        node_tissue_gtex.tsv is absent (e.g. when re-using old TSVs).
        """
        if (self.data_dir / self._TISSUE_FILE).exists():
            tissue_df = self._read(self._TISSUE_FILE)
            required = {"tissue_id", "tissue_name"}
            if required.issubset(set(tissue_df.columns)):
                for r in tissue_df.itertuples(index=False):
                    yield (r.tissue_id, "tissue", {
                        "name": _sanitize(r.tissue_name),
                    })
                return

        # Fallback: derive from edge file
        edge_df = self._read(self._EDGE_FILE)
        for r in (
            edge_df[["tissue_id", "tissue_name"]]
            .drop_duplicates(subset=["tissue_id"])
            .itertuples(index=False)
        ):
            yield (r.tissue_id, "tissue", {
                "name": _sanitize(r.tissue_name),
            })

    # ── Edges ──────────────────────────────────────────────────────── #

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        df = self._read(self._EDGE_FILE).dropna(subset=["gene_symbol", "tissue_id"])

        has_tpm = "median_tpm" in df.columns

        for r in df.itertuples(index=False):
            props: dict = {}
            if has_tpm:
                props["median_tpm"] = float(r.median_tpm)
            gene_id = f"{self._ID_PREFIX}{r.gene_symbol}"
            yield (gene_id, r.tissue_id, self._EDGE_LABEL, props)
