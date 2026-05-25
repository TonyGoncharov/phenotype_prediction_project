"""
src/knowledge_graph/adapters/gene_ontology_adapter.py - adapters for GO annotation edges.

HumanGOAdapter and MouseGOAdapter each read a gene–GO edge TSV produced by the GO
export pipeline and yield gene nodes, GO-term nodes, and gene-has-go-term edges.
"""

from __future__ import annotations

from collections.abc import Generator

import pandas as pd
from src.knowledge_graph.adapters.base import BaseAdapter, EdgeTuple, NodeTuple
from src.knowledge_graph.gene_ids import (
    HUMAN_ORGANISM,
    MOUSE_ORGANISM,
    validate_required_columns,
)

_ASPECT_TO_NAMESPACE: dict[str, str] = {
    'P': 'biological_process',
    'F': 'molecular_function',
    'C': 'cellular_component',
}


class GOAdapter(BaseAdapter):
    """
    Abstract base adapter for GO annotation edges.

    Must not be instantiated directly - use HumanGOAdapter or MouseGOAdapter.
    """

    _GENE_GO_FILE: str = ''
    _GENE_LABEL: str = ''
    _EDGE_LABEL: str = ''
    _ID_PREFIX: str = ''
    _ORGANISM: str = ''
    REQUIRED_COLUMNS: dict[str, set[str]] = {}
    REQUIRED_FILES: list[str] = []

    def __init__(self, data_dir):
        """
        Initialise the adapter and validate required input files.

        Args:
            data_dir: Directory containing processed GO edge TSVs.

        Raises:
            TypeError: If instantiated directly instead of via a subclass.
            FileNotFoundError: If any required TSV file is absent.
        """
        super().__init__(data_dir)
        if not self._GENE_GO_FILE:
            raise TypeError(
                f'{type(self).__name__} must not be instantiated directly - '
                'use HumanGOAdapter or MouseGOAdapter.'
            )
        missing = [f for f in self.REQUIRED_FILES if not (self.data_dir / f).exists()]
        if missing:
            raise FileNotFoundError(
                f'GO layer: missing TSV files in {self.data_dir}:\n  '
                + '\n  '.join(missing)
            )
        self.logger.debug('GO files ready: %s', ', '.join(self.REQUIRED_FILES))
        self._validate_required_columns()

    def _validate_required_columns(self) -> None:
        """
        Validate that all required columns are present in each required file.
        """
        for filename, required in self.REQUIRED_COLUMNS.items():
            df = self._read(filename)
            validate_required_columns(df, required, filename)

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield gene and GO-term nodes for the GO layer.
        """
        yield from self._count_nodes(
            self._raw_nodes(),
            label=self.layer_name,
        )

    def _raw_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield gene nodes then GO-term nodes.
        """
        yield from self._gene_nodes()
        yield from self._go_term_nodes()

    def _gene_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield gene nodes from the gene–GO edge file.
        """
        yield from self._gene_nodes_from_tsv(self._GENE_GO_FILE)

    def _go_term_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield GO-term nodes derived from unique GO IDs in the edge file.
        """
        go_df = self._read(self._GENE_GO_FILE)
        validate_required_columns(go_df, {'go_id', 'aspect'}, self._GENE_GO_FILE)
        aspect_df = go_df[['go_id', 'aspect']].drop_duplicates().dropna()
        go_id_to_ns: dict[str, str] = {
            r.go_id: _ASPECT_TO_NAMESPACE.get(r.aspect, r.aspect)
            for r in aspect_df.itertuples(index=False)
        }
        go_id_to_name: dict[str, str] = {}
        if 'go_name' in go_df.columns:
            name_df = (
                go_df[['go_id', 'go_name']].dropna().drop_duplicates(subset=['go_id'])
            )
            go_id_to_name = dict(zip(name_df['go_id'], name_df['go_name'], strict=True))
        unique_terms = go_df['go_id'].dropna().unique()
        self.logger.debug('GO term nodes to emit: %d', len(unique_terms))
        for go_id in unique_terms:
            yield (
                go_id,
                'go term',
                {
                    'go_id': go_id,
                    'namespace': go_id_to_ns.get(go_id, 'unknown'),
                    'name': self._sanitize(go_id_to_name.get(go_id, '')),
                },
            )

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield gene-has-go-term edges.
        """
        yield from self._count_edges(
            self._raw_edges(),
            label=self.layer_name,
        )

    def _raw_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield gene→GO edges with aggregated evidence codes and aspects.
        """
        df = self._read(self._GENE_GO_FILE).dropna(subset=['gene_id', 'go_id'])
        has_evidence = 'evidence_code' in df.columns
        has_aspect = 'aspect' in df.columns

        agg: dict = {}
        group_cols = ['gene_id', 'go_id']
        if 'gene_symbol' in df.columns:
            group_cols.insert(1, 'gene_symbol')
        if has_evidence:
            agg['evidence_code'] = lambda x: sorted(x.dropna().unique().tolist())
        if has_aspect:
            agg['aspect'] = 'first'

        grouped = (
            df.groupby(group_cols).agg(agg).reset_index()
            if agg
            else df[group_cols].drop_duplicates()
        )
        self.logger.debug('GO edges to emit: %d', len(grouped))

        for r in grouped.itertuples(index=False):
            props: dict = {}
            if has_evidence and (codes := r.evidence_code):
                props['evidence_codes'] = codes
            if has_aspect and pd.notna(getattr(r, 'aspect', None)):
                props['aspect'] = r.aspect
            yield (r.gene_id, r.go_id, self._EDGE_LABEL, props)


class HumanGOAdapter(GOAdapter):
    """
    BioCypher adapter for human gene–GO annotation edges.
    """

    layer_name = 'go'
    _GENE_GO_FILE = 'edge_human_gene_has_go.tsv'
    _GENE_LABEL = 'human gene'
    _EDGE_LABEL = 'human gene has go term'
    _ORGANISM = HUMAN_ORGANISM
    REQUIRED_COLUMNS = {_GENE_GO_FILE: {'gene_id', 'gene_symbol', 'go_id', 'aspect'}}
    REQUIRED_FILES = [_GENE_GO_FILE]


class MouseGOAdapter(GOAdapter):
    """
    BioCypher adapter for mouse gene–GO annotation edges.
    """

    layer_name = 'go'
    _GENE_GO_FILE = 'edge_mouse_gene_has_go.tsv'
    _GENE_LABEL = 'mouse gene'
    _EDGE_LABEL = 'mouse gene has go term'
    _ORGANISM = MOUSE_ORGANISM
    REQUIRED_COLUMNS = {_GENE_GO_FILE: {'gene_id', 'gene_symbol', 'go_id', 'aspect'}}
    REQUIRED_FILES = [_GENE_GO_FILE]
