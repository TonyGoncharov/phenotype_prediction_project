"""
src/knowledge_graph/adapters/gene_expression_adapter.py - adapters for gene expression edges.

HumanExpressionAdapter reads GTEx-derived TSVs; MouseExpressionAdapter reads GXD-derived
TSVs. Both yield gene nodes, tissue nodes, and gene-expressed-in-tissue edges.
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


class HumanExpressionAdapter(BaseAdapter):
    """
    BioCypher adapter for human gene expression edges from GTEx.
    """

    layer_name = 'expression'

    _EDGE_FILE = 'edge_human_gene_expressed_in_tissue.tsv'
    _TISSUE_FILE = 'node_tissue_gtex.tsv'
    _GENE_LABEL = 'human gene'
    _EDGE_LABEL = 'human gene expressed in tissue'
    _ORGANISM = HUMAN_ORGANISM

    REQUIRED_FILES: list[str] = [_EDGE_FILE]
    REQUIRED_COLUMNS: dict[str, set[str]] = {
        _EDGE_FILE: {'gene_id', 'gene_symbol', 'tissue_id', 'tissue_name', 'median_tpm'}
    }

    def __init__(self, data_dir):
        """
        Initialise the adapter and validate required input files.

        Args:
            data_dir: Directory containing processed expression TSVs.

        Raises:
            FileNotFoundError: If any required TSV file is absent.
        """
        super().__init__(data_dir)
        missing = [f for f in self.REQUIRED_FILES if not (self.data_dir / f).exists()]
        if missing:
            raise FileNotFoundError(
                f'Expression layer: missing TSV files in {self.data_dir}:\n  '
                + '\n  '.join(missing)
            )
        tissue_present = (self.data_dir / self._TISSUE_FILE).exists()
        self.logger.debug(
            'Initialised - edge file: %s, tissue node file: %s',
            self._EDGE_FILE,
            'found' if tissue_present else 'absent (will derive from edge file)',
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
        Yield gene and tissue nodes for the expression layer.
        """
        yield from self._count_nodes(
            self._raw_nodes(),
            label=self.layer_name,
        )

    def _raw_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield gene nodes then tissue nodes.
        """
        yield from self._gene_nodes()
        yield from self._tissue_nodes()

    def _gene_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield gene nodes from the edge file.
        """
        yield from self._gene_nodes_from_tsv(self._EDGE_FILE)

    def _tissue_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield tissue nodes from the dedicated tissue file, or derive from the edge file.
        """
        if (self.data_dir / self._TISSUE_FILE).exists():
            tissue_df = self._read(self._TISSUE_FILE)
            required = {'tissue_id', 'tissue_name'}
            validate_required_columns(tissue_df, required, self._TISSUE_FILE)
            self.logger.debug(
                'Tissue nodes to emit: %d (from %s)',
                len(tissue_df),
                self._TISSUE_FILE,
            )
            for r in tissue_df.itertuples(index=False):
                yield (
                    r.tissue_id,
                    'tissue',
                    {
                        'name': self._sanitize(r.tissue_name),
                    },
                )
            return

        self.logger.debug(
            'Tissue node file absent - deriving tissue nodes from edge file'
        )
        edge_df = self._read(self._EDGE_FILE)
        tissue_df = edge_df[['tissue_id', 'tissue_name']].drop_duplicates(
            subset=['tissue_id']
        )
        self.logger.debug('Tissue nodes to emit: %d (derived)', len(tissue_df))
        for r in tissue_df.itertuples(index=False):
            yield (
                r.tissue_id,
                'tissue',
                {
                    'name': self._sanitize(r.tissue_name),
                },
            )

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield gene-expressed-in-tissue edges.
        """
        yield from self._count_edges(
            self._raw_edges(),
            label=self.layer_name,
        )

    def _raw_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield raw expression edges from the edge file.
        """
        df = self._read(self._EDGE_FILE).dropna(subset=['gene_id', 'tissue_id'])
        self.logger.debug('Expression edges to emit: %d', len(df))
        has_tpm = 'median_tpm' in df.columns
        for r in df.itertuples(index=False):
            props: dict = {}
            if has_tpm:
                props['median_tpm'] = float(r.median_tpm)
            yield (r.gene_id, r.tissue_id, self._EDGE_LABEL, props)


class MouseExpressionAdapter(HumanExpressionAdapter):
    """
    BioCypher adapter for mouse gene expression edges from GXD RNA-seq data.
    """

    layer_name = 'expression'

    _EDGE_FILE = 'edge_mouse_gene_expressed_in_tissue.tsv'
    _TISSUE_FILE = 'node_tissue_gxd.tsv'
    _GENE_LABEL = 'mouse gene'
    _EDGE_LABEL = 'mouse gene expressed in tissue'
    _ORGANISM = MOUSE_ORGANISM

    REQUIRED_FILES = [_EDGE_FILE]
    REQUIRED_COLUMNS = {
        _EDGE_FILE: {
            'gene_id',
            'gene_symbol',
            'mgi_gene_id',
            'tissue_id',
            'tissue_name',
            'median_tpm',
            'tpm_level',
            'n_experiments',
        }
    }

    def _tissue_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield mouse tissue nodes, deriving from edge file if dedicated file is absent.
        """
        if (self.data_dir / self._TISSUE_FILE).exists():
            tissue_df = self._read(self._TISSUE_FILE)
            required = {'tissue_id', 'tissue_name'}
            validate_required_columns(tissue_df, required, self._TISSUE_FILE)
        else:
            self.logger.debug(
                'Tissue node file absent - deriving tissue nodes from edge file'
            )
            edge_df = self._read(self._EDGE_FILE)
            tissue_df = edge_df[['tissue_id', 'tissue_name']].drop_duplicates(
                subset=['tissue_id']
            )

        self.logger.debug('Mouse tissue nodes to emit: %d', len(tissue_df))
        for r in tissue_df.itertuples(index=False):
            yield (
                r.tissue_id,
                'tissue',
                {
                    'uberon_id': r.tissue_id,
                    'name': self._sanitize(r.tissue_name),
                },
            )

    def _raw_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield raw mouse expression edges with TPM level and experiment count.
        """
        df = self._read(self._EDGE_FILE).dropna(subset=['gene_id', 'tissue_id'])
        self.logger.debug('Mouse expression edges to emit: %d', len(df))
        for r in df.itertuples(index=False):
            props: dict = {
                'median_tpm': float(r.median_tpm),
                'tpm_level': r.tpm_level,
            }
            if pd.notna(r.n_experiments) and str(r.n_experiments).strip():
                props['n_experiments'] = int(r.n_experiments)
            yield (r.gene_id, r.tissue_id, self._EDGE_LABEL, props)
