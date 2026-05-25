"""
src/knowledge_graph/adapters/gene_to_phenotype_adapter.py - adapters for gene→MP phenotype edges.

HumanPhenotypeAdapter and MousePhenotypeAdapter read gene–MP top-term edge TSVs
produced by the phenotype export pipeline and yield gene nodes, MP top-term nodes,
and gene-has-mp-top-term edges.
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


class PhenotypeAdapter(BaseAdapter):
    """
    Abstract base adapter for gene→MP phenotype edges.

    Must not be instantiated directly - use HumanPhenotypeAdapter or MousePhenotypeAdapter.
    """

    _EDGE_FILE: str = ''
    _GENE_LABEL: str = ''
    _EDGE_LABEL: str = ''
    _ORGANISM: str = ''
    REQUIRED_COLUMNS: dict[str, set[str]] = {}
    REQUIRED_FILES: list[str] = []

    def __init__(self, data_dir):
        """
        Initialise the adapter and validate required input files.

        Args:
            data_dir: Directory containing processed phenotype edge TSVs.

        Raises:
            TypeError: If instantiated directly instead of via a subclass.
            FileNotFoundError: If any required TSV file is absent.
        """
        super().__init__(data_dir)
        if not self._EDGE_FILE:
            raise TypeError(
                f'{type(self).__name__} must not be instantiated directly - '
                'use HumanPhenotypeAdapter or MousePhenotypeAdapter.'
            )
        missing = [f for f in self.REQUIRED_FILES if not (self.data_dir / f).exists()]
        if missing:
            raise FileNotFoundError(
                f'Phenotype layer: missing TSV files in {self.data_dir}:\n  '
                + '\n  '.join(missing)
            )
        self.logger.debug(
            'Initialised - edge file: %s, MP names: %s',
            self._EDGE_FILE,
            'found' if (self.data_dir / 'node_mp_top_names.tsv').exists() else 'absent',
        )
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
        Yield gene and MP top-term nodes for the phenotype layer.
        """
        yield from self._count_nodes(
            self._unique_nodes(self._all_nodes()),
            label=self.layer_name,
        )

    def _all_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield gene nodes then MP top-term nodes.
        """
        yield from self._gene_nodes()
        yield from self._mp_top_term_nodes()

    def _gene_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield gene nodes from the edge file.
        """
        yield from self._gene_nodes_from_tsv(self._EDGE_FILE)

    def _mp_top_term_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield MP top-term nodes derived from unique mp_top_id values in the edge file.
        """
        df = self._read(self._EDGE_FILE)
        validate_required_columns(df, {'mp_top_id'}, self._EDGE_FILE)
        mp_names = self._load_mp_names()
        unique_tops = df['mp_top_id'].dropna().unique()
        self.logger.debug('MP top-term nodes to emit: %d', len(unique_tops))
        for mp_top_id in unique_tops:
            yield (
                mp_top_id,
                'mp top term',
                {
                    'mp_id': mp_top_id,
                    'name': self._sanitize(mp_names.get(mp_top_id, '')),
                },
            )

    def _load_mp_names(self) -> dict[str, str]:
        """
        Load MP top-term names from node_mp_top_names.tsv if present.

        Returns:
            Dict mapping mp_id → name, or empty dict if the file is absent.
        """
        names_path = self.data_dir / 'node_mp_top_names.tsv'
        if not names_path.exists():
            return {}
        names_df = self._read('node_mp_top_names.tsv')
        if {'mp_id', 'name'}.issubset(names_df.columns):
            return dict(zip(names_df['mp_id'], names_df['name'], strict=True))
        return {}

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield gene-has-mp-top-term edges.
        """
        yield from self._count_edges(self._raw_edges(), label=self.layer_name)

    def _raw_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield raw phenotype edges with optional HP source IDs and distance.
        """
        df = self._read(self._EDGE_FILE).dropna(subset=['gene_id', 'mp_top_id'])

        has_hp_ids = 'source_hp_ids' in df.columns
        has_dist = 'min_hp_to_anchor_dist' in df.columns

        for r in df.itertuples(index=False):
            props: dict = {'derived': True}
            if has_hp_ids and pd.notna(getattr(r, 'source_hp_ids', None)):
                props['source_hp_ids'] = str(r.source_hp_ids).split('|')
            if has_dist and pd.notna(getattr(r, 'min_hp_to_anchor_dist', None)):
                props['min_hp_to_anchor_distance'] = int(r.min_hp_to_anchor_dist)
            yield (r.gene_id, r.mp_top_id, self._EDGE_LABEL, props)


class HumanPhenotypeAdapter(PhenotypeAdapter):
    """
    BioCypher adapter for human gene→MP top-term phenotype edges.
    """

    layer_name = 'phenotype'
    _EDGE_FILE = 'edge_human_gene_has_mp_top.tsv'
    _GENE_LABEL = 'human gene'
    _EDGE_LABEL = 'human gene has mp top term'
    _ORGANISM = HUMAN_ORGANISM
    REQUIRED_COLUMNS = {_EDGE_FILE: {'gene_id', 'gene_symbol', 'mp_top_id'}}
    REQUIRED_FILES = [_EDGE_FILE]


class MousePhenotypeAdapter(PhenotypeAdapter):
    """
    BioCypher adapter for mouse gene→MP top-term phenotype edges.
    """

    layer_name = 'phenotype'
    _EDGE_FILE = 'edge_mouse_gene_has_mp_top.tsv'
    _GENE_LABEL = 'mouse gene'
    _EDGE_LABEL = 'mouse gene has mp top term'
    _ORGANISM = MOUSE_ORGANISM
    REQUIRED_COLUMNS = {_EDGE_FILE: {'gene_id', 'gene_symbol', 'mp_top_id'}}
    REQUIRED_FILES = [_EDGE_FILE]
