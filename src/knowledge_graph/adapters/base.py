"""
src/knowledge_graph/adapters/base.py - shared base class for all BioCypher adapters.

Provides TSV reading, node/edge deduplication, count logging, and the
`_gene_nodes_from_tsv` helper that all per-layer adapters reuse.
"""

from __future__ import annotations

import logging
from collections.abc import Generator, Iterator
from pathlib import Path

import pandas as pd
from src.knowledge_graph.gene_ids import validate_required_columns

NodeTuple = tuple[str, str, dict]
EdgeTuple = tuple[str, str, str, dict]


class BaseAdapter:
    """
    Abstract base class for all BioCypher knowledge-graph adapters.

    Subclasses implement get_nodes() and get_edges() to yield tuples for
    BioCypher. Helper methods handle TSV reading, deduplication, and logging.
    """

    layer_name: str = ''

    def __init__(self, data_dir: str | Path):
        """
        Initialise the adapter with the directory containing processed TSVs.

        Args:
            data_dir: Path to the directory holding per-layer TSV files.
        """
        self.data_dir = Path(data_dir)
        self.logger = logging.getLogger(
            f'{type(self).__module__}.{type(self).__name__}'
        )

    def _read(self, filename: str) -> pd.DataFrame:
        """
        Read a tab-separated TSV from data_dir and return a DataFrame.

        Args:
            filename: Filename relative to data_dir.

        Raises:
            ValueError: If the file is completely empty (no rows or header).
        """
        path = self.data_dir / filename
        self.logger.debug('Reading TSV: %s', path)
        try:
            df = pd.read_csv(path, sep='\t', dtype=str)
            self.logger.debug('Loaded %d rows from %s', len(df), filename)
            return df
        except pd.errors.EmptyDataError:
            raise ValueError(
                f'TSV file is completely empty (no content at all): {path}'
            ) from None

    @staticmethod
    def _unique_nodes(
        tuples: Iterator[NodeTuple],
    ) -> Generator[NodeTuple, None, None]:
        """
        Yield unique nodes by node_id, discarding duplicates.

        Args:
            tuples: Iterator of (node_id, label, props) tuples.

        Yields:
            NodeTuple for each unique node_id encountered.
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
        """
        Yield unique edges by (src, tgt, label) key, discarding duplicates.

        Args:
            tuples: Iterator of (src, tgt, label, props) tuples.

        Yields:
            EdgeTuple for each unique (src, tgt, label) combination.
        """
        seen: set[tuple[str, str, str]] = set()
        for src, tgt, label, props in tuples:
            key = (src, tgt, label)
            if key not in seen:
                seen.add(key)
                yield src, tgt, label, props

    def _count_nodes(
        self,
        gen: Iterator[NodeTuple],
        label: str = '',
    ) -> Generator[NodeTuple, None, None]:
        """
        Pass nodes through while logging counts by label.

        Args:
            gen: Source node iterator.
            label: Optional prefix for log messages.

        Yields:
            All NodeTuples from *gen* unchanged.
        """
        counts: dict[str, int] = {}
        for node_id, node_label, props in gen:
            counts[node_label] = counts.get(node_label, 0) + 1
            yield node_id, node_label, props
        total = sum(counts.values())
        tag = f'[{label}] ' if label else ''
        if counts:
            breakdown = ', '.join(f'{lbl}: {n}' for lbl, n in sorted(counts.items()))
            self.logger.info('%sNodes emitted - total: %d  (%s)', tag, total, breakdown)
        else:
            self.logger.warning('%sNo nodes emitted.', tag)

    def _count_edges(
        self,
        gen: Iterator[EdgeTuple],
        label: str = '',
    ) -> Generator[EdgeTuple, None, None]:
        """
        Pass edges through while logging counts by edge label.

        Args:
            gen: Source edge iterator.
            label: Optional prefix for log messages.

        Yields:
            All EdgeTuples from *gen* unchanged.
        """
        counts: dict[str, int] = {}
        for src, tgt, edge_label, props in gen:
            counts[edge_label] = counts.get(edge_label, 0) + 1
            yield src, tgt, edge_label, props
        total = sum(counts.values())
        tag = f'[{label}] ' if label else ''
        if counts:
            breakdown = ', '.join(f'{lbl}: {n}' for lbl, n in sorted(counts.items()))
            self.logger.info('%sEdges emitted - total: %d  (%s)', tag, total, breakdown)
        else:
            self.logger.warning('%sNo edges emitted.', tag)

    @staticmethod
    def _sanitize(value: str) -> str:
        """
        Strip single-quote characters from a string value.

        Args:
            value: Input string.

        Returns:
            String with all single-quote characters removed.
        """
        return value.replace("'", '')

    def _check_empty_df(self, df: pd.DataFrame, name: str) -> bool:
        """
        Log a warning and return True if *df* is empty.

        Args:
            df: DataFrame to check.
            name: Human-readable name for log messages.

        Returns:
            True if df is empty, False otherwise.
        """
        if df.empty:
            self.logger.warning('%s file is empty.', name)
            return True
        return False

    def _gene_nodes_from_tsv(
        self,
        tsv_filename: str,
        gene_col: str = 'gene_id',
    ) -> Generator[NodeTuple, None, None]:
        """
        Yield gene node tuples from a TSV file.

        Args:
            tsv_filename: Filename relative to data_dir.
            gene_col: Column containing the gene ID.

        Yields:
            NodeTuple for each unique gene_id in the TSV.
        """
        df = self._read(tsv_filename)
        validate_required_columns(df, {gene_col}, f'{tsv_filename} gene nodes')
        unique_genes = df[gene_col].dropna().unique()
        self.logger.debug('Gene nodes to emit: %d', len(unique_genes))
        gene_label: str = getattr(self, '_GENE_LABEL', 'gene')
        organism: str = getattr(self, '_ORGANISM', '')
        symbol_by_id: dict[str, str] = {}
        if 'gene_symbol' in df.columns:
            symbol_by_id = (
                df[[gene_col, 'gene_symbol']]
                .dropna()
                .drop_duplicates(subset=[gene_col])
                .set_index(gene_col)['gene_symbol']
                .to_dict()
            )
        for gene_id in unique_genes:
            props = {'organism': organism}
            if symbol := symbol_by_id.get(gene_id):
                props['symbol'] = symbol
            yield (gene_id, gene_label, props)

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield all node tuples for this adapter layer.

        Yields:
            NodeTuple objects consumed by BioCypher.
        """
        yield from ()

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield all edge tuples for this adapter layer.

        Yields:
            EdgeTuple objects consumed by BioCypher.
        """
        yield from ()
