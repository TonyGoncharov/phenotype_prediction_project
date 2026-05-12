"""uberon_adapter.py - Uberon anatomical ontology adapter for the gene-phenotype KG.

Reads the two TSV files produced by uberon_export.py and yields:
  Nodes : uberon term — one per Uberon anatomical term referenced by a GTEx tissue
  Edges : tissue mapped to uberon — GTEx tissue → Uberon term
"""

from __future__ import annotations

from typing import Generator

from src.adapters.base import BaseAdapter, NodeTuple, EdgeTuple


class UberonAdapter(BaseAdapter):
    """Read Uberon anatomical ontology nodes and tissue-mapping edges.

    Do not instantiate directly — use HumanUberonAdapter or MouseUberonAdapter.

    Expected files in data_dir
    --------------------------
    node_uberon_terms.tsv
        uberon_id | name | definition | synonyms

    edge_tissue_mapped_to_uberon.tsv
        gtex_tissue_id | uberon_id | gtex_tissue_name
    """

    _UBERON_NODE_FILE: str = "node_uberon_terms.tsv"
    _UBERON_EDGE_FILE: str = "edge_tissue_mapped_to_uberon.tsv"

    def __init__(self, data_dir):
        super().__init__(data_dir)
        if type(self) is UberonAdapter:
            raise TypeError(
                "UberonAdapter must not be instantiated directly — "
                "use HumanUberonAdapter or MouseUberonAdapter."
            )

        self._has_nodes = (self.data_dir / self._UBERON_NODE_FILE).exists()
        self._has_edges = (self.data_dir / self._UBERON_EDGE_FILE).exists()

        if not self._has_nodes:
            self.logger.warning(
                "%s not found — Uberon term nodes will be omitted. "
                "Run the Uberon export step first.",
                self._UBERON_NODE_FILE,
            )
        if not self._has_edges:
            self.logger.warning(
                "%s not found — tissue→Uberon edges will be omitted. "
                "Run the Uberon export step first.",
                self._UBERON_EDGE_FILE,
            )

    # ── Nodes ──────────────────────────────────────────────────────────── #

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        if not self._has_nodes:
            return
        yield from self._count_nodes(
            self._uberon_term_nodes(),
            label=self.layer_name,
        )

    def _uberon_term_nodes(self) -> Generator[NodeTuple, None, None]:
        df = self._read(self._UBERON_NODE_FILE)
        if self._check_empty_df(df, "Uberon node"):
            return

        self.logger.debug("Uberon term nodes to emit: %d", len(df))
        for row in df.itertuples(index=False):
            uberon_id: str = row.uberon_id
            if not uberon_id or not uberon_id.startswith("UBERON:"):
                self.logger.debug("Skipping non-Uberon ID: %s", uberon_id)
                continue

            props: dict = {
                "uberon_id": uberon_id,
                "name":      getattr(row, "name", "") or "",
            }
            # Optional columns written by the export step
            if hasattr(row, "definition") and row.definition:
                props["definition"] = row.definition
            if hasattr(row, "synonyms") and row.synonyms:
                # synonyms are stored pipe-separated in the TSV
                props["synonyms"] = [s.strip() for s in str(row.synonyms).split("|") if s.strip()]

            yield (uberon_id, "uberon term", props)

    # ── Edges ──────────────────────────────────────────────────────────── #

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        if not self._has_edges:
            return
        yield from self._count_edges(
            self._tissue_to_uberon_edges(),
            label=self.layer_name,
        )

    def _tissue_to_uberon_edges(self) -> Generator[EdgeTuple, None, None]:
        df = self._read(self._UBERON_EDGE_FILE)
        if self._check_empty_df(df, "Uberon edge"):
            return

        required = {"gtex_tissue_id", "uberon_id"}
        missing_cols = required - set(df.columns)
        if missing_cols:
            self.logger.error(
                "Uberon edge TSV missing required columns: %s", missing_cols
            )
            return

        df = df.dropna(subset=["gtex_tissue_id", "uberon_id"])
        self.logger.debug("Tissue→Uberon edges to emit: %d", len(df))

        for row in df.itertuples(index=False):
            gtex_id: str  = row.gtex_tissue_id
            uberon_id: str = row.uberon_id

            if not uberon_id.startswith("UBERON:"):
                self.logger.debug("Skipping non-Uberon target: %s", uberon_id)
                continue

            props: dict = {}
            if hasattr(row, "gtex_tissue_name") and row.gtex_tissue_name:
                props["gtex_tissue_name"] = row.gtex_tissue_name

            yield (gtex_id, uberon_id, "tissue mapped to uberon", props)


# ── Species-specific subclasses ───────────────────────────────────────────── #

class HumanUberonAdapter(UberonAdapter):
    """Uberon adapter for the human graph; reads the same shared TSV files as MouseUberonAdapter."""
    layer_name = "uberon"


class MouseUberonAdapter(UberonAdapter):
    """Uberon anatomical ontology adapter for the mouse graph.

    Note: the mouse expression layer is not yet implemented, so tissue nodes
    for the mouse graph are not currently present.  This adapter is provided
    for forward-compatibility; it will produce empty output until mouse tissue
    data is added to the graph.
    """
    layer_name = "uberon"
