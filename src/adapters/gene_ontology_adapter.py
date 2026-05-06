"""src/adapters/gene_ontology_adapter.py"""

from __future__ import annotations

from typing import Generator

import pandas as pd

from src.adapters.base import BaseAdapter, NodeTuple, EdgeTuple

_ASPECT_TO_NAMESPACE: dict[str, str] = {
    "P": "biological_process",
    "F": "molecular_function",
    "C": "cellular_component",
}


class GOAdapter(BaseAdapter):
    """Read GO annotation edges for one species.

    Do not instantiate directly — use HumanGOAdapter or MouseGOAdapter.
    """

    _GENE_GO_FILE: str = ""
    _GENE_LABEL: str = ""
    _EDGE_LABEL: str = ""
    _ID_PREFIX: str = ""
    _ORGANISM: str = ""

    def __init__(self, data_dir):
        super().__init__(data_dir)
        if not self._GENE_GO_FILE:
            raise TypeError(
                f"{type(self).__name__} must not be instantiated directly — "
                "use HumanGOAdapter or MouseGOAdapter."
            )
        if not (self.data_dir / self._GENE_GO_FILE).exists():
            raise FileNotFoundError(
                f"GO layer: {self._GENE_GO_FILE} not found in {self.data_dir}.\n"
                "Run the GO export step first."
            )
        self.logger.debug("GO annotation file found: %s", self._GENE_GO_FILE)

    # ── Nodes ──────────────────────────────────────────────────────── #

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        yield from self._count_nodes(
            self._raw_nodes(),
            label=self.layer_name,
        )

    def _raw_nodes(self) -> Generator[NodeTuple, None, None]:
        yield from self._gene_nodes()
        yield from self._go_term_nodes()

    def _gene_nodes(self) -> Generator[NodeTuple, None, None]:
        df = self._read(self._GENE_GO_FILE)
        unique_genes = df["gene_symbol"].dropna().unique()
        self.logger.debug("Gene nodes to emit: %d", len(unique_genes))
        for symbol in unique_genes:
            node_id = f"{self._ID_PREFIX}{symbol}" if self._ID_PREFIX else symbol
            yield (node_id, self._GENE_LABEL, {
                "symbol": symbol,
                "organism": self._ORGANISM,
            })

    def _go_term_nodes(self) -> Generator[NodeTuple, None, None]:
        """Yield one GoTerm node per unique GO ID."""
        go_df = self._read(self._GENE_GO_FILE)
        aspect_df = go_df[["go_id", "aspect"]].drop_duplicates().dropna()
        go_id_to_ns: dict[str, str] = {
            r.go_id: _ASPECT_TO_NAMESPACE.get(r.aspect, r.aspect)
            for r in aspect_df.itertuples(index=False)
        }
        go_id_to_name: dict[str, str] = {}
        if "go_name" in go_df.columns:
            name_df = go_df[["go_id", "go_name"]].dropna().drop_duplicates(subset=["go_id"])
            go_id_to_name = dict(zip(name_df["go_id"], name_df["go_name"]))
        unique_terms = go_df["go_id"].dropna().unique()
        self.logger.debug("GO term nodes to emit: %d", len(unique_terms))
        for go_id in unique_terms:
            yield (go_id, "go term", {
                "go_id": go_id,
                "namespace": go_id_to_ns.get(go_id, "unknown"),
                "name": self._sanitize(go_id_to_name.get(go_id, "")),
            })

    # ── Edges ──────────────────────────────────────────────────────── #

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        yield from self._count_edges(
            self._raw_edges(),
            label=self.layer_name,
        )

    def _raw_edges(self) -> Generator[EdgeTuple, None, None]:
        """Yield gene → GO term edges, aggregating evidence codes per (gene, term) pair."""
        df = self._read(self._GENE_GO_FILE).dropna(subset=["gene_symbol", "go_id"])
        has_evidence = "evidence_code" in df.columns
        has_aspect   = "aspect" in df.columns

        agg: dict = {}
        if has_evidence:
            agg["evidence_code"] = lambda x: sorted(x.dropna().unique().tolist())
        if has_aspect:
            agg["aspect"] = "first"

        grouped = (
            df.groupby(["gene_symbol", "go_id"]).agg(agg).reset_index()
            if agg
            else df[["gene_symbol", "go_id"]].drop_duplicates()
        )
        self.logger.debug("GO edges to emit: %d", len(grouped))

        for r in grouped.itertuples(index=False):
            props: dict = {}
            if has_evidence and (codes := r.evidence_code):
                props["evidence_codes"] = codes
            if has_aspect and pd.notna(getattr(r, "aspect", None)):
                props["aspect"] = r.aspect
            gene_id = f"{self._ID_PREFIX}{r.gene_symbol}" if self._ID_PREFIX else r.gene_symbol
            yield (gene_id, r.go_id, self._EDGE_LABEL, props)


# ── Species-specific subclasses ────────────────────────────────────── #

class HumanGOAdapter(GOAdapter):
    """GO annotation adapter for the human graph."""
    layer_name    = "go"
    _GENE_GO_FILE = "edge_human_gene_has_go.tsv"
    _GENE_LABEL   = "human gene"
    _EDGE_LABEL   = "human gene has go term"
    _ID_PREFIX    = "HGNC:"
    _ORGANISM     = "Homo sapiens"


class MouseGOAdapter(GOAdapter):
    """GO annotation adapter for the mouse graph."""
    layer_name    = "go"
    _GENE_GO_FILE = "edge_mouse_gene_has_go.tsv"
    _GENE_LABEL   = "mouse gene"
    _EDGE_LABEL   = "mouse gene has go term"
    _ID_PREFIX    = ""
    _ORGANISM     = "Mus musculus"
