"""
src/knowledge_graph/adapters/ppi_adapter.py - adapters for protein-protein interaction edges.

HumanPPIAdapter and MousePPIAdapter each read three TSVs (protein nodes,
gene-encodes-protein edges, protein-interacts-with edges) produced by the PPI
export pipeline and yield the corresponding BioCypher tuples.
"""

from __future__ import annotations

from collections.abc import Generator

from src.knowledge_graph.adapters.base import BaseAdapter, EdgeTuple, NodeTuple
from src.knowledge_graph.gene_ids import validate_required_columns


def _parse_pipe_list(value: object) -> list[str]:
    """
    Split a pipe-delimited string into a list; return empty list for non-strings.

    Args:
        value: Input value (expected to be a string or None/NaN).

    Returns:
        List of non-empty substrings, or [] if value is not a non-empty string.
    """
    if not isinstance(value, str) or not value.strip():
        return []
    return value.split('|')


class PPIAdapter(BaseAdapter):
    """
    Abstract base adapter for protein-protein interaction edges.

    Must not be instantiated directly - use HumanPPIAdapter or MousePPIAdapter.
    """

    _PROTEIN_FILE: str = ''
    _ENCODES_FILE: str = ''
    _PPI_FILE: str = ''
    _ENCODES_LABEL: str = ''
    _PPI_LABEL: str = ''
    REQUIRED_COLUMNS: dict[str, set[str]] = {}
    REQUIRED_FILES: list[str] = []

    def __init__(self, data_dir):
        """
        Initialise the adapter and validate required input files.

        Args:
            data_dir: Directory containing processed PPI TSVs.

        Raises:
            TypeError: If instantiated directly.
            FileNotFoundError: If any required TSV file is absent.
        """
        super().__init__(data_dir)
        if type(self) is PPIAdapter:
            raise TypeError('Use HumanPPIAdapter or MousePPIAdapter.')
        missing = [f for f in self.REQUIRED_FILES if not (self.data_dir / f).exists()]
        if missing:
            raise FileNotFoundError(
                f'PPI layer: missing TSV files in {self.data_dir}:\n  '
                + '\n  '.join(missing)
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
        Yield protein nodes for the PPI layer.
        """
        yield from self._count_nodes(self._protein_nodes(), label=self.layer_name)

    def _protein_nodes(self) -> Generator[NodeTuple, None, None]:
        """
        Yield protein nodes from the protein node file.
        """
        df = self._read(self._PROTEIN_FILE)
        if self._check_empty_df(df, 'Protein node'):
            return
        for row in df.itertuples(index=False):
            yield (
                row.protein_id,
                'protein',
                {
                    'gene_symbol': row.gene_symbol,
                    'organism': row.organism,
                },
            )

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield gene-encodes-protein and protein-interacts-with edges.
        """
        yield from self._count_edges(self._raw_edges(), label=self.layer_name)

    def _raw_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield all edges: gene→protein encodes edges then protein↔protein PPI edges.
        """
        yield from self._gene_encodes_edges()
        yield from self._ppi_edges()

    def _gene_encodes_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield gene-encodes-protein edges from the encodes edge file.
        """
        df = self._read(self._ENCODES_FILE)
        if self._check_empty_df(df, 'Gene-encodes-protein'):
            return
        for row in df.itertuples(index=False):
            yield (row.gene_id, row.protein_id, self._ENCODES_LABEL, {})

    def _ppi_edges(self) -> Generator[EdgeTuple, None, None]:
        """
        Yield protein-interacts-with edges with experimental metadata.
        """
        df = self._read(self._PPI_FILE)
        if self._check_empty_df(df, 'PPI edge'):
            return

        has_systems = 'experimental_systems' in df.columns
        has_sys_types = 'experimental_system_types' in df.columns
        has_throughput = 'throughput' in df.columns
        has_pubmed = 'pubmed_ids' in df.columns
        has_source = 'source_db' in df.columns

        for row in df.itertuples(index=False):
            props: dict = {}
            if has_systems:
                props['experimental_systems'] = _parse_pipe_list(
                    getattr(row, 'experimental_systems', '')
                )
            if has_sys_types:
                props['experimental_system_types'] = _parse_pipe_list(
                    getattr(row, 'experimental_system_types', '')
                )
            if has_throughput and (tp := getattr(row, 'throughput', None)):
                props['throughput'] = tp
            if has_pubmed:
                props['pubmed_ids'] = _parse_pipe_list(getattr(row, 'pubmed_ids', ''))
            if has_source and (src := getattr(row, 'source_db', None)):
                props['source_db'] = src
            yield (row.protein_id_a, row.protein_id_b, self._PPI_LABEL, props)


class HumanPPIAdapter(PPIAdapter):
    """
    BioCypher adapter for human protein-protein interaction edges from BioGRID.
    """

    layer_name = 'ppi'
    _PROTEIN_FILE = 'node_human_protein.tsv'
    _ENCODES_FILE = 'edge_human_gene_encodes_protein.tsv'
    _PPI_FILE = 'edge_human_ppi.tsv'
    _ENCODES_LABEL = 'human gene encodes protein'
    _PPI_LABEL = 'human protein interacts with'
    REQUIRED_FILES = [_PROTEIN_FILE, _ENCODES_FILE, _PPI_FILE]
    REQUIRED_COLUMNS = {
        _PROTEIN_FILE: {'protein_id', 'gene_symbol', 'organism'},
        _ENCODES_FILE: {'gene_id', 'gene_symbol', 'protein_id'},
        _PPI_FILE: {'protein_id_a', 'protein_id_b'},
    }


class MousePPIAdapter(PPIAdapter):
    """
    BioCypher adapter for mouse protein-protein interaction edges from BioGRID.
    """

    layer_name = 'ppi'
    _PROTEIN_FILE = 'node_mouse_protein.tsv'
    _ENCODES_FILE = 'edge_mouse_gene_encodes_protein.tsv'
    _PPI_FILE = 'edge_mouse_ppi.tsv'
    _ENCODES_LABEL = 'mouse gene encodes protein'
    _PPI_LABEL = 'mouse protein interacts with'
    REQUIRED_FILES = [_PROTEIN_FILE, _ENCODES_FILE, _PPI_FILE]
    REQUIRED_COLUMNS = {
        _PROTEIN_FILE: {'protein_id', 'gene_symbol', 'organism'},
        _ENCODES_FILE: {'gene_id', 'gene_symbol', 'protein_id'},
        _PPI_FILE: {'protein_id_a', 'protein_id_b'},
    }
