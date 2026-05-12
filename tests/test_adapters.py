"""Tests for all BioCypher adapter classes.

Covers:
  - TypeError when instantiating abstract base adapters directly
  - FileNotFoundError when required TSV files are missing
  - Successful init and node/edge emission with minimal fixture data
"""

import pytest

from src.adapters.gene_to_phenotype_adapter import (
    HumanPhenotypeAdapter,
    MousePhenotypeAdapter,
    PhenotypeAdapter,
)
from src.adapters.gene_expression_adapter import HumanExpressionAdapter
from src.adapters.gene_ontology_adapter import (
    GOAdapter,
    HumanGOAdapter,
    MouseGOAdapter,
)
from src.adapters.ppi_adapter import HumanPPIAdapter
from src.adapters.uberon_adapter import (
    HumanUberonAdapter,
    MouseUberonAdapter,
    UberonAdapter,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────

def _write_tsv(path, header: str, rows: list[str]) -> None:
    path.write_text("\n".join([header] + rows) + "\n")


@pytest.fixture()
def human_pheno_dir(tmp_path):
    _write_tsv(
        tmp_path / "edge_human_gene_has_mp_top.tsv",
        "gene_symbol\tmp_top_id\tsource_hp_ids\tmin_hp_to_anchor_dist",
        ["TP53\tMP:0005369\tHP:0000001\t1", "BRCA1\tMP:0005369\tHP:0000002\t2"],
    )
    return tmp_path


@pytest.fixture()
def mouse_pheno_dir(tmp_path):
    _write_tsv(
        tmp_path / "edge_mouse_gene_has_mp_top.tsv",
        "gene_symbol\tmp_top_id",
        ["Trp53\tMP:0005369", "Brca1\tMP:0005369"],
    )
    return tmp_path


@pytest.fixture()
def human_go_dir(tmp_path):
    _write_tsv(
        tmp_path / "edge_human_gene_has_go.tsv",
        "gene_symbol\tgo_id\taspect\tevidence_code\tgo_name",
        ["TP53\tGO:0008150\tP\tIEA\tapoptotic process"],
    )
    return tmp_path


@pytest.fixture()
def mouse_go_dir(tmp_path):
    _write_tsv(
        tmp_path / "edge_mouse_gene_has_go.tsv",
        "gene_symbol\tgo_id\taspect\tevidence_code",
        ["Trp53\tGO:0008150\tP\tIEA"],
    )
    return tmp_path


@pytest.fixture()
def expression_dir(tmp_path):
    _write_tsv(
        tmp_path / "edge_human_gene_expressed_in_tissue.tsv",
        "gene_symbol\ttissue_id\ttissue_name\tmedian_tpm",
        ["TP53\tGTEX:Lung\tLung\t12.5"],
    )
    return tmp_path


@pytest.fixture()
def ppi_dir(tmp_path):
    _write_tsv(
        tmp_path / "node_human_protein.tsv",
        "protein_id\tgene_symbol\torganism",
        ["PROTEIN:TP53\tTP53\tHomo sapiens"],
    )
    _write_tsv(
        tmp_path / "edge_human_gene_encodes_protein.tsv",
        "gene_symbol\tprotein_id",
        ["TP53\tPROTEIN:TP53"],
    )
    _write_tsv(
        tmp_path / "edge_human_ppi.tsv",
        "protein_id_a\tprotein_id_b",
        ["PROTEIN:TP53\tPROTEIN:BRCA1"],
    )
    return tmp_path


@pytest.fixture()
def uberon_dir(tmp_path):
    _write_tsv(
        tmp_path / "node_uberon_terms.tsv",
        "uberon_id\tname",
        ["UBERON:0002048\tlung"],
    )
    _write_tsv(
        tmp_path / "edge_tissue_mapped_to_uberon.tsv",
        "gtex_tissue_id\tuberon_id\tgtex_tissue_name",
        ["GTEX:Lung\tUBERON:0002048\tLung"],
    )
    return tmp_path


# ── PhenotypeAdapter ──────────────────────────────────────────────────────────

class TestPhenotypeAdapterBase:
    def test_raises_typeerror_when_instantiated_directly(self, tmp_path):
        with pytest.raises(TypeError, match="must not be instantiated directly"):
            PhenotypeAdapter(tmp_path)


class TestHumanPhenotypeAdapter:
    def test_raises_on_missing_edge_file(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Phenotype layer"):
            HumanPhenotypeAdapter(tmp_path)

    def test_init_succeeds_with_valid_files(self, human_pheno_dir):
        adapter = HumanPhenotypeAdapter(human_pheno_dir)
        assert adapter.layer_name == "phenotype"

    def test_required_files_attribute(self):
        assert "edge_human_gene_has_mp_top.tsv" in HumanPhenotypeAdapter.REQUIRED_FILES

    def test_get_nodes_emits_gene_and_mp_nodes(self, human_pheno_dir):
        adapter = HumanPhenotypeAdapter(human_pheno_dir)
        nodes = list(adapter.get_nodes())
        labels = {label for _, label, _ in nodes}
        assert "human gene" in labels
        assert "mp top term" in labels

    def test_gene_nodes_have_hgnc_prefix(self, human_pheno_dir):
        adapter = HumanPhenotypeAdapter(human_pheno_dir)
        gene_ids = [nid for nid, label, _ in adapter.get_nodes() if label == "human gene"]
        assert all(gid.startswith("HGNC:") for gid in gene_ids)

    def test_get_edges_emits_edges(self, human_pheno_dir):
        adapter = HumanPhenotypeAdapter(human_pheno_dir)
        edges = list(adapter.get_edges())
        assert len(edges) > 0
        assert all(label == "human gene has mp top term" for _, _, label, _ in edges)


class TestMousePhenotypeAdapter:
    def test_raises_on_missing_edge_file(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Phenotype layer"):
            MousePhenotypeAdapter(tmp_path)

    def test_init_succeeds_with_valid_files(self, mouse_pheno_dir):
        adapter = MousePhenotypeAdapter(mouse_pheno_dir)
        assert adapter.layer_name == "phenotype"

    def test_required_files_attribute(self):
        assert "edge_mouse_gene_has_mp_top.tsv" in MousePhenotypeAdapter.REQUIRED_FILES

    def test_gene_nodes_have_no_prefix(self, mouse_pheno_dir):
        adapter = MousePhenotypeAdapter(mouse_pheno_dir)
        gene_ids = [nid for nid, label, _ in adapter.get_nodes() if label == "mouse gene"]
        assert all(not gid.startswith("HGNC:") for gid in gene_ids)

    def test_get_edges_emits_edges(self, mouse_pheno_dir):
        adapter = MousePhenotypeAdapter(mouse_pheno_dir)
        edges = list(adapter.get_edges())
        assert len(edges) > 0


# ── GOAdapter ─────────────────────────────────────────────────────────────────

class TestGOAdapterBase:
    def test_raises_typeerror_when_instantiated_directly(self, tmp_path):
        with pytest.raises(TypeError, match="must not be instantiated directly"):
            GOAdapter(tmp_path)


class TestHumanGOAdapter:
    def test_raises_on_missing_go_file(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="GO layer"):
            HumanGOAdapter(tmp_path)

    def test_init_succeeds_with_valid_files(self, human_go_dir):
        adapter = HumanGOAdapter(human_go_dir)
        assert adapter.layer_name == "go"

    def test_required_files_attribute(self):
        assert "edge_human_gene_has_go.tsv" in HumanGOAdapter.REQUIRED_FILES

    def test_get_nodes_emits_gene_and_go_nodes(self, human_go_dir):
        adapter = HumanGOAdapter(human_go_dir)
        nodes = list(adapter.get_nodes())
        labels = {label for _, label, _ in nodes}
        assert "human gene" in labels
        assert "go term" in labels

    def test_gene_nodes_have_hgnc_prefix(self, human_go_dir):
        adapter = HumanGOAdapter(human_go_dir)
        gene_ids = [nid for nid, label, _ in adapter.get_nodes() if label == "human gene"]
        assert all(gid.startswith("HGNC:") for gid in gene_ids)

    def test_get_edges_emit_go_edges(self, human_go_dir):
        adapter = HumanGOAdapter(human_go_dir)
        edges = list(adapter.get_edges())
        assert len(edges) > 0
        assert all(label == "human gene has go term" for _, _, label, _ in edges)


class TestMouseGOAdapter:
    def test_raises_on_missing_go_file(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="GO layer"):
            MouseGOAdapter(tmp_path)

    def test_init_succeeds_with_valid_files(self, mouse_go_dir):
        adapter = MouseGOAdapter(mouse_go_dir)
        assert adapter.layer_name == "go"

    def test_required_files_attribute(self):
        assert "edge_mouse_gene_has_go.tsv" in MouseGOAdapter.REQUIRED_FILES

    def test_gene_nodes_have_no_prefix(self, mouse_go_dir):
        adapter = MouseGOAdapter(mouse_go_dir)
        gene_ids = [nid for nid, label, _ in adapter.get_nodes() if label == "mouse gene"]
        assert all(not gid.startswith("HGNC:") for gid in gene_ids)


# ── HumanExpressionAdapter ───────────────────────────────────────────────────

class TestHumanExpressionAdapter:
    def test_raises_on_missing_edge_file(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Expression layer"):
            HumanExpressionAdapter(tmp_path)

    def test_init_succeeds_with_valid_files(self, expression_dir):
        adapter = HumanExpressionAdapter(expression_dir)
        assert adapter.layer_name == "expression"

    def test_required_files_attribute(self):
        assert "edge_human_gene_expressed_in_tissue.tsv" in HumanExpressionAdapter.REQUIRED_FILES

    def test_get_nodes_emits_gene_and_tissue_nodes(self, expression_dir):
        adapter = HumanExpressionAdapter(expression_dir)
        nodes = list(adapter.get_nodes())
        labels = {label for _, label, _ in nodes}
        assert "human gene" in labels
        assert "tissue" in labels

    def test_gene_nodes_have_hgnc_prefix(self, expression_dir):
        adapter = HumanExpressionAdapter(expression_dir)
        gene_ids = [nid for nid, label, _ in adapter.get_nodes() if label == "human gene"]
        assert all(gid.startswith("HGNC:") for gid in gene_ids)

    def test_get_edges_emit_expression_edges(self, expression_dir):
        adapter = HumanExpressionAdapter(expression_dir)
        edges = list(adapter.get_edges())
        assert len(edges) > 0
        assert all(label == "human gene expressed in tissue" for _, _, label, _ in edges)

    def test_edge_properties_contain_median_tpm(self, expression_dir):
        adapter = HumanExpressionAdapter(expression_dir)
        edges = list(adapter.get_edges())
        assert "median_tpm" in edges[0][3]


# ── HumanPPIAdapter ───────────────────────────────────────────────────────────

class TestHumanPPIAdapter:
    def test_raises_on_missing_all_files(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="PPI layer"):
            HumanPPIAdapter(tmp_path)

    def test_raises_on_partial_missing_files(self, tmp_path):
        _write_tsv(
            tmp_path / "node_human_protein.tsv",
            "protein_id\tgene_symbol\torganism",
            ["PROTEIN:TP53\tTP53\tHomo sapiens"],
        )
        with pytest.raises(FileNotFoundError, match="PPI layer"):
            HumanPPIAdapter(tmp_path)

    def test_init_succeeds_with_all_files(self, ppi_dir):
        adapter = HumanPPIAdapter(ppi_dir)
        assert adapter.layer_name == "ppi"

    def test_required_files_attribute_is_public(self):
        assert hasattr(HumanPPIAdapter, "REQUIRED_FILES")
        assert not hasattr(HumanPPIAdapter, "_REQUIRED_FILES")

    def test_required_files_lists_all_three_tsvs(self):
        assert len(HumanPPIAdapter.REQUIRED_FILES) == 3

    def test_get_nodes_emits_protein_nodes(self, ppi_dir):
        adapter = HumanPPIAdapter(ppi_dir)
        nodes = list(adapter.get_nodes())
        labels = {label for _, label, _ in nodes}
        assert "protein" in labels

    def test_get_edges_emits_both_edge_types(self, ppi_dir):
        adapter = HumanPPIAdapter(ppi_dir)
        edges = list(adapter.get_edges())
        edge_labels = {label for _, _, label, _ in edges}
        assert "human gene encodes protein" in edge_labels
        assert "human protein interacts with" in edge_labels


# ── UberonAdapter ─────────────────────────────────────────────────────────────

class TestUberonAdapterBase:
    def test_raises_typeerror_when_instantiated_directly(self, tmp_path):
        with pytest.raises(TypeError, match="must not be instantiated directly"):
            UberonAdapter(tmp_path)


class TestHumanUberonAdapter:
    def test_init_succeeds_without_files(self, tmp_path):
        """Uberon uses silent degradation — missing files are not an error."""
        adapter = HumanUberonAdapter(tmp_path)
        assert adapter.layer_name == "uberon"

    def test_empty_output_when_files_missing(self, tmp_path):
        adapter = HumanUberonAdapter(tmp_path)
        assert list(adapter.get_nodes()) == []
        assert list(adapter.get_edges()) == []

    def test_get_nodes_emits_uberon_nodes_when_file_present(self, uberon_dir):
        adapter = HumanUberonAdapter(uberon_dir)
        nodes = list(adapter.get_nodes())
        assert any(label == "uberon term" for _, label, _ in nodes)

    def test_get_edges_emits_mapping_edges_when_file_present(self, uberon_dir):
        adapter = HumanUberonAdapter(uberon_dir)
        edges = list(adapter.get_edges())
        assert any(label == "tissue mapped to uberon" for _, _, label, _ in edges)


class TestMouseUberonAdapter:
    def test_init_succeeds_without_files(self, tmp_path):
        adapter = MouseUberonAdapter(tmp_path)
        assert adapter.layer_name == "uberon"

    def test_empty_output_when_files_missing(self, tmp_path):
        adapter = MouseUberonAdapter(tmp_path)
        assert list(adapter.get_nodes()) == []
        assert list(adapter.get_edges()) == []
