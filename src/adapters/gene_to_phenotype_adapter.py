'''
src/adapters/gene_to_phenotype_adapter.py
'''

from __future__ import annotations

from typing import Generator

from src.adapters.base import BaseAdapter, NodeTuple, EdgeTuple


def _sanitize(value: str) -> str:
    """Remove single quotes that break neo4j-admin CSV import (quote char = ')."""
    return value.replace("'", "")


class PhenotypeAdapter(BaseAdapter):
    _EDGE_FILE: str = ""
    _GENE_LABEL: str = ""
    _EDGE_LABEL: str = ""
    _ID_PREFIX: str = ""
    _ORGANISM: str = ""

    REQUIRED_FILES: list[str] = []

    def __init__(self, data_dir):
        super().__init__(data_dir)
        if not self._EDGE_FILE:
            raise TypeError(
                f"{type(self).__name__} must not be instantiated directly — "
                "use HumanPhenotypeAdapter or MousePhenotypeAdapter."
            )
        missing = [f for f in self.REQUIRED_FILES if not (self.data_dir / f).exists()]
        if missing:
            raise FileNotFoundError(
                f"Phenotype layer: missing TSV files in {self.data_dir}:\n  "
                + "\n  ".join(missing)
            )

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        yield from self._gene_nodes()
        yield from self._mp_top_term_nodes()

    def _gene_nodes(self) -> Generator[NodeTuple, None, None]:
        df = self._read(self._EDGE_FILE)
        for symbol in df["gene_symbol"].dropna().unique():
            node_id = f"{self._ID_PREFIX}{symbol}" if self._ID_PREFIX else symbol
            yield (node_id, self._GENE_LABEL, {
                "symbol": symbol,
                "organism": self._ORGANISM,
            })

    def _mp_top_term_nodes(self) -> Generator[NodeTuple, None, None]:
        df = self._read(self._EDGE_FILE)
        # Load human-readable names produced by the phenotype export step
        names_file = "node_mp_top_names.tsv"
        mp_names: dict[str, str] = {}
        if (self.data_dir / names_file).exists():
            names_df = self._read(names_file)
            if "mp_id" in names_df.columns and "name" in names_df.columns:
                mp_names = dict(zip(names_df["mp_id"], names_df["name"]))
        for mp_top_id in df["mp_top_id"].dropna().unique():
            yield (mp_top_id, "mp top term", {
                "mp_id": mp_top_id,
                "name": _sanitize(mp_names.get(mp_top_id, "")),
            })

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        df = self._read(self._EDGE_FILE).dropna(subset=["gene_symbol", "mp_top_id"])
        for r in df.itertuples(index=False):
            gene_id = f"{self._ID_PREFIX}{r.gene_symbol}" if self._ID_PREFIX else r.gene_symbol
            yield (gene_id, r.mp_top_id, self._EDGE_LABEL, {})


class HumanPhenotypeAdapter(PhenotypeAdapter):
    layer_name  = "phenotype"
    _EDGE_FILE  = "edge_human_gene_has_mp_top.tsv"
    _GENE_LABEL = "human gene"
    _EDGE_LABEL = "human gene has mp top term"
    _ID_PREFIX  = "HGNC:"
    _ORGANISM   = "Homo sapiens"

    REQUIRED_FILES = [_EDGE_FILE]


class MousePhenotypeAdapter(PhenotypeAdapter):
    layer_name  = "phenotype"
    _EDGE_FILE  = "edge_mouse_gene_has_mp_top.tsv"
    _GENE_LABEL = "mouse gene"
    _EDGE_LABEL = "mouse gene has mp top term"
    _ID_PREFIX  = ""
    _ORGANISM   = "Mus musculus"

    REQUIRED_FILES = [_EDGE_FILE]