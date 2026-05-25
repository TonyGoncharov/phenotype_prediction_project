"""
src/knowledge_graph/adapters/uberon_adapter.py - adapter for Uberon anatomical ontology edges.

HumanUberonAdapter reads the tissue-to-Uberon mapping and Uberon term node files
produced by the Uberon export pipeline and yields tissue-mapped-to-uberon edges
and uberon-term nodes.
"""

from __future__ import annotations

from collections.abc import Generator

from src.knowledge_graph.adapters.base import BaseAdapter, EdgeTuple, NodeTuple
from src.knowledge_graph.gene_ids import validate_required_columns


class UberonAdapter(BaseAdapter):
    """
    Abstract base adapter for Uberon anatomical ontology edges.

    Must not be instantiated directly - use HumanUberonAdapter.
    """

    _UBERON_NODE_FILE: str = 'node_uberon_terms.tsv'
    _UBERON_EDGE_FILE: str = 'edge_tissue_mapped_to_uberon.tsv'

    def __init__(self, data_dir):
        """
        Initialise the adapter and check for at least one Uberon file.

        Args:
            data_dir: Directory containing processed Uberon TSVs.

        Raises:
            TypeError: If instantiated directly.
            FileNotFoundError: If neither Uberon node nor edge file is present.
        """
        super().__init__(data_dir)
        if type(self) is UberonAdapter:
            raise TypeError(
                'UberonAdapter must not be instantiated directly - '
                'use HumanUberonAdapter.'
            )

        self._has_nodes = (self.data_dir / self._UBERON_NODE_FILE).exists()
        self._has_edges = (self.data_dir / self._UBERON_EDGE_FILE).exists()

        if not self._has_nodes and not self._has_edges:
            raise FileNotFoundError(
                f'Uberon layer: missing TSV files in {self.data_dir}:\n  '
                + '\n  '.join([self._UBERON_NODE_FILE, self._UBERON_EDGE_FILE])
            )

        if not self._has_nodes:
            self.logger.warning(
                '%s not found - Uberon term nodes will be omitted. '
                'Run the Uberon export step first.',
                self._UBERON_NODE_FILE,
            )
        if not self._has_edges:
            self.logger.warning(
                '%s not found - tissue→Uberon edges will be omitted. '
                'Run the Uberon export step first.',
                self._UBERON_EDGE_FILE,
            )

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield Uberon term nodes if the node file is present.
        """
        if not self._has_nodes:
            return
        yield from self._count_nodes(
            self._uberon_term_nodes(),
            label=self.layer_name,
        )

    def _uberon_term_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield Uberon term nodes from the node file.
        """
        df = self._read(self._UBERON_NODE_FILE)
        if self._check_empty_df(df, 'Uberon node'):
            return
        validate_required_columns(df, {'uberon_id'}, self._UBERON_NODE_FILE)

        self.logger.debug('Uberon term nodes to emit: %d', len(df))
        for row in df.itertuples(index=False):
            uberon_id: str = row.uberon_id
            if not uberon_id or not uberon_id.startswith('UBERON:'):
                self.logger.debug('Skipping non-Uberon ID: %s', uberon_id)
                continue

            props: dict = {
                'uberon_id': uberon_id,
                'name': getattr(row, 'name', '') or '',
            }

            if hasattr(row, 'definition') and row.definition:
                props['definition'] = row.definition
            if hasattr(row, 'synonyms') and row.synonyms:
                props['synonyms'] = [
                    s.strip() for s in str(row.synonyms).split('|') if s.strip()
                ]

            yield (uberon_id, 'uberon term', props)

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield tissue-mapped-to-uberon edges if the edge file is present.
        """
        if not self._has_edges:
            return
        yield from self._count_edges(
            self._tissue_to_uberon_edges(),
            label=self.layer_name,
        )

    def _tissue_to_uberon_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield tissue→Uberon mapping edges from the edge file.
        """
        df = self._read(self._UBERON_EDGE_FILE)
        if self._check_empty_df(df, 'Uberon edge'):
            return

        validate_required_columns(
            df, {'gtex_tissue_id', 'uberon_id'}, self._UBERON_EDGE_FILE
        )

        df = df.dropna(subset=['gtex_tissue_id', 'uberon_id'])
        self.logger.debug('Tissue→Uberon edges to emit: %d', len(df))

        for row in df.itertuples(index=False):
            gtex_id: str = row.gtex_tissue_id
            uberon_id: str = row.uberon_id

            if not uberon_id.startswith('UBERON:'):
                self.logger.debug('Skipping non-Uberon target: %s', uberon_id)
                continue

            props: dict = {}
            if hasattr(row, 'gtex_tissue_name') and row.gtex_tissue_name:
                props['gtex_tissue_name'] = row.gtex_tissue_name

            yield (gtex_id, uberon_id, 'tissue mapped to uberon', props)


class HumanUberonAdapter(UberonAdapter):
    """
    BioCypher adapter for the human Uberon anatomical ontology layer.
    """

    layer_name = 'uberon'


class MouseUberonAdapter(UberonAdapter):
    """
    Placeholder - mouse expression TSVs embed UBERON IDs directly; no separate layer.
    """

    layer_name = 'uberon'

    def __init__(self, data_dir):
        """
        Always raises NotImplementedError.

        Raises:
            NotImplementedError: Mouse Uberon mapping is not a separate build layer.
        """
        raise NotImplementedError(
            'Mouse Uberon mapping is not a separate build layer. '
            'Mouse expression TSVs use UBERON tissue IDs directly.'
        )
