"""
src/adapters/gene_ontology_adapter.py
"""

from __future__ import annotations

import itertools
from typing import Generator

import pandas as pd

from src.adapters.base import BaseAdapter, NodeTuple, EdgeTuple

_ASPECT_TO_NAMESPACE: dict[str, str] = {
    "P": "biological_process",
    "F": "molecular_function",
    "C": "cellular_component",
}


def _sanitize(value: str) -> str:
    """Remove single quotes that break neo4j-admin CSV import (quote char = ')."""
    return value.replace("'", "")


class GOAdapter(BaseAdapter):
    """Read GO annotation edges for one species.

    Do not instantiate directly — use HumanGOAdapter or MouseGOAdapter.
    """

    # ── Subclasses configure these ─────────────────────────────────── #

    #: e.g. "edge_human_gene_has_go.tsv"
    _GENE_GO_FILE: str = ""

    #: BioCypher node label for genes, e.g. "human gene"
    _GENE_LABEL: str = ""

    #: BioCypher edge label, e.g. "human gene has go term"
    _EDGE_LABEL: str = ""

    #: Prepended to gene_symbol to form the gene node ID.
    _ID_PREFIX: str = ""

    #: Human-readable species name stored as a node property.
    _ORGANISM: str = ""

    def __init__(self, data_dir):
        super().__init__(data_dir)
        if not self._GENE_GO_FILE:
            raise TypeError(
                f"{type(self).__name__} must not be instantiated directly — "
                "use HumanGOAdapter or MouseGOAdapter."
            )
        self._has_go = (self.data_dir / self._GENE_GO_FILE).exists()
        if not self._has_go:
            print(
                f"  [{type(self).__name__}] NOTE: {self._GENE_GO_FILE} not found — "
                "GO nodes/edges will be omitted.  Run the GO export step first."
            )

    # ── Nodes ──────────────────────────────────────────────────────── #

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        if not self._has_go:
            return
        yield from self._gene_nodes()
        yield from self._go_term_nodes()

    def _gene_nodes(self) -> Generator[NodeTuple, None, None]:
        df = self._read(self._GENE_GO_FILE)
        for symbol in df["gene_symbol"].dropna().unique():
            node_id = f"{self._ID_PREFIX}{symbol}" if self._ID_PREFIX else symbol
            yield (node_id, self._GENE_LABEL, {
                "symbol": symbol,
                "organism": self._ORGANISM,
            })

    def _go_term_nodes(self) -> Generator[NodeTuple, None, None]:
        """Yield one GoTerm node per unique GO ID.

        namespace is read from the aspect column (P/F/C → full name).
        name is read from the go_name column when available.
        """
        go_df = self._read(self._GENE_GO_FILE)
        aspect_df = go_df[["go_id", "aspect"]].drop_duplicates().dropna()
        go_id_to_ns: dict[str, str] = {
            r.go_id: _ASPECT_TO_NAMESPACE.get(r.aspect, r.aspect)
            for r in aspect_df.itertuples(index=False)
        }
        # Build go_id → human-readable name map (first non-null wins)
        go_id_to_name: dict[str, str] = {}
        if "go_name" in go_df.columns:
            name_df = go_df[["go_id", "go_name"]].dropna().drop_duplicates(subset=["go_id"])
            go_id_to_name = dict(zip(name_df["go_id"], name_df["go_name"]))
        for go_id in go_df["go_id"].dropna().unique():
            yield (go_id, "go term", {
                "go_id": go_id,
                "namespace": go_id_to_ns.get(go_id, "unknown"),
                "name": _sanitize(go_id_to_name.get(go_id, "")),
            })

    # ── Edges ──────────────────────────────────────────────────────── #

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        if not self._has_go:
            return
        yield from self._gene_has_go_edges()

    def _gene_has_go_edges(self) -> Generator[EdgeTuple, None, None]:
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