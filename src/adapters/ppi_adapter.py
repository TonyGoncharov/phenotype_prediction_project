"""ppi_adapter.py – BioCypher adapter for human protein-protein interactions.

Reads three TSV files produced by ppi_export.py and emits:

Nodes
-----
* ``protein``               — one node per human protein
                              ID: ``PROTEIN:<HGNC_symbol>``

Edges
-----
* ``human gene encodes protein``   — Gene → Protein (gene-to-protein mapping)
* ``human protein interacts with`` — Protein ↔ Protein (physical PPI)

Layer name: ``ppi``
"""

from __future__ import annotations

from typing import Generator

from src.adapters.base import BaseAdapter, EdgeTuple, NodeTuple

# File names written by ppi_export.py
_PROTEIN_NODES_FILE        = "node_human_protein.tsv"
_GENE_ENCODES_FILE         = "edge_human_gene_encodes_protein.tsv"
_PPI_EDGES_FILE            = "edge_human_ppi.tsv"

_GENE_ID_PREFIX            = "HGNC:"   # must match HumanPhenotypeAdapter / HumanGOAdapter
_PROTEIN_ID_PREFIX         = "PROTEIN:"


def _parse_pipe_list(value: str) -> list[str]:
    """Convert a pipe-separated string back to a list (e.g. ``'a|b|c'`` → ``['a','b','c']``)."""
    if not isinstance(value, str) or not value.strip():
        return []
    return value.split("|")


class HumanPPIAdapter(BaseAdapter):
    """Adapter for human protein-protein interaction data (BioGRID).

    Expects three TSV files in *data_dir* (written by ``ppi_export.py``):

    * ``node_human_protein.tsv``
    * ``edge_human_gene_encodes_protein.tsv``
    * ``edge_human_ppi.tsv``
    """

    layer_name = "ppi"
    _REQUIRED_FILES: list[str] = [_PROTEIN_NODES_FILE, _GENE_ENCODES_FILE, _PPI_EDGES_FILE]

    def __init__(self, data_dir):
        super().__init__(data_dir)
        missing = [f for f in self._REQUIRED_FILES if not (self.data_dir / f).exists()]
        if missing:
            raise FileNotFoundError(
                f"PPI layer: missing TSV files in {self.data_dir}:\n  "
                + "\n  ".join(missing)
            )
        self.logger.debug("PPI files ready: %s", ", ".join(self._REQUIRED_FILES))

    # ── Nodes ──────────────────────────────────────────────────────────

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        yield from self._count_nodes(self._protein_nodes(), label=self.layer_name)

    def _protein_nodes(self) -> Generator[NodeTuple, None, None]:
        df = self._read(_PROTEIN_NODES_FILE)
        if df.empty:
            self.logger.warning("Protein node file is empty.")
            return
        self.logger.debug("Protein nodes to emit: %d", len(df))
        for row in df.itertuples(index=False):
            yield (
                row.protein_id,
                "protein",
                {
                    "gene_symbol": row.gene_symbol,
                    "organism":    row.organism,
                },
            )

    # ── Edges ──────────────────────────────────────────────────────────

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        yield from self._count_edges(self._raw_edges(), label=self.layer_name)

    def _raw_edges(self) -> Generator[EdgeTuple, None, None]:
        yield from self._gene_encodes_edges()
        yield from self._ppi_edges()

    def _gene_encodes_edges(self) -> Generator[EdgeTuple, None, None]:
        """Yield Gene → Protein 'encodes' edges."""
        df = self._read(_GENE_ENCODES_FILE)
        if df.empty:
            self.logger.warning("Gene-encodes-protein file is empty.")
            return
        self.logger.debug("Gene→Protein edges to emit: %d", len(df))
        for row in df.itertuples(index=False):
            gene_id = f"{_GENE_ID_PREFIX}{row.gene_symbol}"
            yield (
                gene_id,
                row.protein_id,
                "human gene encodes protein",
                {},
            )

    def _ppi_edges(self) -> Generator[EdgeTuple, None, None]:
        """Yield Protein ↔ Protein physical interaction edges."""
        df = self._read(_PPI_EDGES_FILE)
        if df.empty:
            self.logger.warning("PPI edge file is empty.")
            return
        self.logger.debug("PPI edges to emit: %d", len(df))

        has_systems    = "experimental_systems"       in df.columns
        has_sys_types  = "experimental_system_types"  in df.columns
        has_throughput = "throughput"                 in df.columns
        has_pubmed     = "pubmed_ids"                 in df.columns
        has_source     = "source_db"                  in df.columns

        for row in df.itertuples(index=False):
            props: dict = {}
            if has_systems:
                props["experimental_systems"] = _parse_pipe_list(
                    getattr(row, "experimental_systems", "")
                )
            if has_sys_types:
                props["experimental_system_types"] = _parse_pipe_list(
                    getattr(row, "experimental_system_types", "")
                )
            if has_throughput and (tp := getattr(row, "throughput", None)):
                props["throughput"] = tp
            if has_pubmed:
                props["pubmed_ids"] = _parse_pipe_list(
                    getattr(row, "pubmed_ids", "")
                )
            if has_source and (src := getattr(row, "source_db", None)):
                props["source_db"] = src

            yield (
                row.protein_id_a,
                row.protein_id_b,
                "human protein interacts with",
                props,
            )
