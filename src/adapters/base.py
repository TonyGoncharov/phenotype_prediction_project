"""
src/adapters/base.py - shared base class for all layer adapters.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Generator, Iterator

import pandas as pd

# Concrete type aliases instead of bare `tuple`
NodeTuple = tuple[str, str, dict]        # (id, label, properties)
EdgeTuple = tuple[str, str, str, dict]   # (source_id, target_id, label, properties)


class BaseAdapter:
    """
    Base class for all layer adapters.

    Provides:
      - _read()            — read a TSV from data_dir into a DataFrame
      - _unique_nodes()    — deduplicate node tuples by id
      - _count_nodes()     — wrap a node generator and log the final count
      - _count_edges()     — wrap an edge generator and log the final count

    Each subclass must:
      - Set ``layer_name`` (used by build_graph.py for --skip-layers)
      - Implement get_nodes() and/or get_edges()
      - Accept __init__(self, data_dir) with no other required arguments
        so that build_graph.py can instantiate all adapters uniformly.
    """

    layer_name: str = ""

    def __init__(self, data_dir: str | Path):
        self.data_dir = Path(data_dir)
        self.logger = logging.getLogger(
            f"{type(self).__module__}.{type(self).__name__}"
        )

    # ── I/O ────────────────────────────────────────────────────────── #

    def _read(self, filename: str) -> pd.DataFrame:
        """Read a tab-separated file from data_dir.  Returns empty DataFrame on empty file."""
        path = self.data_dir / filename
        self.logger.debug("Reading TSV: %s", path)
        try:
            df = pd.read_csv(path, sep="\t", dtype=str)
            self.logger.debug("Loaded %d rows from %s", len(df), filename)
            return df
        except pd.errors.EmptyDataError:
            raise ValueError(f"TSV file is completely empty (no content at all): {path}") from None

    # ── Deduplication ──────────────────────────────────────────────── #

    @staticmethod
    def _unique_nodes(
        tuples: Iterator[NodeTuple],
    ) -> Generator[NodeTuple, None, None]:
        """
        Yield node tuples, skipping duplicates by id.

        Always pass a lazy iterator (e.g. itertools.chain), not a materialised
        tuple/list — unpacking generators with (*gen1, *gen2) loads everything
        into memory and defeats the purpose of lazy generation.
        """
        seen: set[str] = set()
        for node_id, label, props in tuples:
            if node_id not in seen:
                seen.add(node_id)
                yield node_id, label, props

    @staticmethod
    def _unique_edges(
        tuples: Iterator[EdgeTuple],
    ) -> Generator[EdgeTuple, None, None]:
        """Yield edge tuples, skipping duplicates by (source_id, target_id, label)."""
        seen: set[tuple[str, str, str]] = set()
        for src, tgt, label, props in tuples:
            key = (src, tgt, label)
            if key not in seen:
                seen.add(key)
                yield src, tgt, label, props

    # ── Instrumented generators ────────────────────────────────────── #

    def _count_nodes(
        self,
        gen: Iterator[NodeTuple],
        label: str = "",
    ) -> Generator[NodeTuple, None, None]:
        """
        Wrap a node generator, logging a count summary when it is exhausted.
        """
        counts: dict[str, int] = {}
        for node_id, node_label, props in gen:
            counts[node_label] = counts.get(node_label, 0) + 1
            yield node_id, node_label, props
        total = sum(counts.values())
        tag = f"[{label}] " if label else ""
        if counts:
            breakdown = ", ".join(f"{lbl}: {n}" for lbl, n in sorted(counts.items()))
            self.logger.info("%sNodes emitted — total: %d  (%s)", tag, total, breakdown)
        else:
            self.logger.warning("%sNo nodes emitted.", tag)

    def _count_edges(
        self,
        gen: Iterator[EdgeTuple],
        label: str = "",
    ) -> Generator[EdgeTuple, None, None]:
        """
        Wrap an edge generator, logging a count summary when it is exhausted.
        """
        counts: dict[str, int] = {}
        for src, tgt, edge_label, props in gen:
            counts[edge_label] = counts.get(edge_label, 0) + 1
            yield src, tgt, edge_label, props
        total = sum(counts.values())
        tag = f"[{label}] " if label else ""
        if counts:
            breakdown = ", ".join(f"{lbl}: {n}" for lbl, n in sorted(counts.items()))
            self.logger.info("%sEdges emitted — total: %d  (%s)", tag, total, breakdown)
        else:
            self.logger.warning("%sNo edges emitted.", tag)

    # ── String helpers ────────────────────────────────────────────────── #

    @staticmethod
    def _sanitize(value: str) -> str:
        """
        Remove single quotes that break neo4j-admin CSV import (quote char = ').
        """
        return value.replace("'", "")

    # ── DataFrame guards ───────────────────────────────────────────── #

    def _check_empty_df(self, df: "pd.DataFrame", name: str) -> bool:
        """Return True and log a warning if *df* has no rows; False otherwise."""
        if df.empty:
            self.logger.warning("%s file is empty.", name)
            return True
        return False

    # ── Shared node helpers ────────────────────────────────────────── #

    def _gene_nodes_from_tsv(
        self,
        tsv_filename: str,
        gene_col: str = "gene_symbol",
    ) -> Generator[NodeTuple, None, None]:
        """Yield one gene node per unique value in *gene_col* of *tsv_filename*.

        Uses self._ID_PREFIX, self._GENE_LABEL, and self._ORGANISM — all of
        which must be defined by the calling subclass.
        """
        df = self._read(tsv_filename)
        unique_genes = df[gene_col].dropna().unique()
        self.logger.debug("Gene nodes to emit: %d", len(unique_genes))
        id_prefix: str = getattr(self, "_ID_PREFIX", "")
        gene_label: str = getattr(self, "_GENE_LABEL", "gene")
        organism: str = getattr(self, "_ORGANISM", "")
        for symbol in unique_genes:
            node_id = f"{id_prefix}{symbol}" if id_prefix else symbol
            yield (node_id, gene_label, {"symbol": symbol, "organism": organism})

    # ── Default no-op generators ───────────────────────────────────── #

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        yield from ()

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        yield from ()
