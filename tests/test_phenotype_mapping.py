"""Tests for src/layers/phenotype_mapping.py

Covers lift_hp_to_anchors() and build_hp_to_mp_top() without loading
any real OBO files — the mock ontology is a plain dict-backed stub that
satisfies the subset of the pronto API these functions actually call.
"""

import pytest

from src.layers.phenotype_mapping import build_hp_to_mp_top, lift_hp_to_anchors


# ── Minimal pronto-compatible ontology stub ───────────────────────────────── #

class _MockTerm:
    def __init__(self, term_id: str, ancestor_ids: list[str]) -> None:
        self.id = term_id
        self._ancestor_ids = ancestor_ids

    def superclasses(self, with_self: bool = True):
        terms = [_MockTerm(aid, []) for aid in self._ancestor_ids]
        if with_self:
            terms.append(self)
        return terms


class _MockOntology:
    """dict-backed stub for pronto.Ontology.

    Constructor argument: {term_id: [ancestor_id, ...]}
    """
    def __init__(self, terms: dict[str, list[str]]) -> None:
        self._terms = terms

    def __contains__(self, term_id: str) -> bool:
        return term_id in self._terms

    def __getitem__(self, term_id: str) -> _MockTerm:
        return _MockTerm(term_id, self._terms[term_id])


# Shared ontology fixture
#
#   HP:ROOT
#   ├── HP:ANCHOR_A
#   │     └── HP:LEAF      (maps to both ANCHOR_A and ANCHOR_B)
#   ├── HP:ANCHOR_B
#   │     └── HP:LEAF
#   └── HP:ORPHAN           (no anchor ancestor)
#
ONT = _MockOntology({
    "HP:ROOT":     [],
    "HP:ANCHOR_A": ["HP:ROOT"],
    "HP:ANCHOR_B": ["HP:ROOT"],
    "HP:LEAF":     ["HP:ANCHOR_A", "HP:ANCHOR_B", "HP:ROOT"],
    "HP:ORPHAN":   ["HP:ROOT"],
})

ANCHORS = {"HP:ANCHOR_A", "HP:ANCHOR_B"}


# ── lift_hp_to_anchors ────────────────────────────────────────────────────── #

class TestLiftHpToAnchors:
    def test_unknown_term_returns_empty(self):
        assert lift_hp_to_anchors("HP:UNKNOWN", ONT, ANCHORS) == []

    def test_term_is_itself_an_anchor(self):
        assert lift_hp_to_anchors("HP:ANCHOR_A", ONT, ANCHORS) == ["HP:ANCHOR_A"]

    def test_leaf_with_one_anchor_ancestor(self):
        ont = _MockOntology({
            "HP:ROOT":   [],
            "HP:ANCHOR": ["HP:ROOT"],
            "HP:LEAF":   ["HP:ANCHOR", "HP:ROOT"],
        })
        assert lift_hp_to_anchors("HP:LEAF", ont, {"HP:ANCHOR"}) == ["HP:ANCHOR"]

    def test_leaf_with_multiple_anchor_ancestors(self):
        result = lift_hp_to_anchors("HP:LEAF", ONT, ANCHORS)
        assert set(result) == ANCHORS

    def test_result_is_sorted(self):
        result = lift_hp_to_anchors("HP:LEAF", ONT, ANCHORS)
        assert result == sorted(result)

    def test_no_anchor_ancestor_returns_empty(self):
        assert lift_hp_to_anchors("HP:ORPHAN", ONT, ANCHORS) == []

    def test_empty_anchor_set_returns_empty(self):
        assert lift_hp_to_anchors("HP:LEAF", ONT, set()) == []


# ── build_hp_to_mp_top ────────────────────────────────────────────────────── #

class TestBuildHpToMpTop:
    def test_leaf_maps_to_both_mp_terms(self):
        hp2mp = {"HP:ANCHOR_A": ["MP:001"], "HP:ANCHOR_B": ["MP:002"]}
        result = build_hp_to_mp_top(["HP:LEAF"], ONT, hp2mp)
        assert set(result["hp_id"]) == {"HP:LEAF"}
        assert set(result["mp_top_id"]) == {"MP:001", "MP:002"}

    def test_unmapped_term_produces_no_rows(self):
        hp2mp = {"HP:ANCHOR_A": ["MP:001"]}
        result = build_hp_to_mp_top(["HP:ORPHAN"], ONT, hp2mp)
        assert result.empty

    def test_empty_input_returns_empty_dataframe(self):
        result = build_hp_to_mp_top([], ONT, {"HP:ANCHOR_A": ["MP:001"]})
        assert result.empty
        assert list(result.columns) == ["hp_id", "mp_top_id"]

    def test_anchor_term_maps_directly(self):
        hp2mp = {"HP:ANCHOR_A": ["MP:001"]}
        result = build_hp_to_mp_top(["HP:ANCHOR_A"], ONT, hp2mp)
        assert len(result) == 1
        assert result.iloc[0]["hp_id"] == "HP:ANCHOR_A"
        assert result.iloc[0]["mp_top_id"] == "MP:001"

    def test_duplicate_inputs_deduplicated(self):
        hp2mp = {"HP:ANCHOR_A": ["MP:001"]}
        ont = _MockOntology({
            "HP:ROOT":   [],
            "HP:ANCHOR_A": ["HP:ROOT"],
            "HP:LEAF":   ["HP:ANCHOR_A", "HP:ROOT"],
        })
        result = build_hp_to_mp_top(["HP:LEAF", "HP:LEAF"], ont, hp2mp)
        assert len(result) == 1

    def test_output_columns(self):
        result = build_hp_to_mp_top(["HP:LEAF"], ONT, {"HP:ANCHOR_A": ["MP:001"]})
        assert list(result.columns) == ["hp_id", "mp_top_id"]
