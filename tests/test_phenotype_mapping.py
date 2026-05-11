"""Tests for src/layers/phenotype_mapping.py

Covers lift_hp_to_anchors() and build_hp_to_mp_top() without loading
any real OBO files — the mock ontology is a plain dict-backed stub that
satisfies the subset of the pronto API these functions actually call.
"""

import pytest

from src.layers.phenotype_mapping import build_hp_to_mp_top, lift_hp_to_anchors


# ── Minimal pronto-compatible ontology stub ───────────────────────────────── #

class _MockTerm:
    def __init__(self, term_id: str, parent_ids: list[str]) -> None:
        self.id = term_id
        self._parent_ids = parent_ids

    def superclasses(self, with_self: bool = True, distance: int = 0):
        # distance=1 is what the BFS calls — return direct parents only.
        # _parent_ids holds direct parents; that is all BFS ever needs.
        terms = [_MockTerm(pid, []) for pid in self._parent_ids]
        if with_self:
            terms.append(self)
        return terms


class _MockOntology:
    """dict-backed stub for pronto.Ontology.

    Constructor argument: {term_id: [direct_parent_id, ...]}
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
#   │     └── HP:LEAF      (distance 1 to ANCHOR_A, distance 1 to ANCHOR_B)
#   ├── HP:ANCHOR_B
#   │     └── HP:LEAF
#   └── HP:ORPHAN           (no anchor ancestor)
#
# _parent_ids stores DIRECT parents only (BFS uses distance=1).
ONT = _MockOntology({
    "HP:ROOT":     [],
    "HP:ANCHOR_A": ["HP:ROOT"],
    "HP:ANCHOR_B": ["HP:ROOT"],
    "HP:LEAF":     ["HP:ANCHOR_A", "HP:ANCHOR_B"],
    "HP:ORPHAN":   ["HP:ROOT"],
})

ANCHORS = {"HP:ANCHOR_A", "HP:ANCHOR_B"}


# ── lift_hp_to_anchors ────────────────────────────────────────────────────── #

class TestLiftHpToAnchors:
    def test_unknown_term_returns_empty(self):
        assert lift_hp_to_anchors("HP:UNKNOWN", ONT, ANCHORS) == {}

    def test_term_is_itself_an_anchor(self):
        assert lift_hp_to_anchors("HP:ANCHOR_A", ONT, ANCHORS) == {"HP:ANCHOR_A": 0}

    def test_leaf_with_one_anchor_ancestor(self):
        ont = _MockOntology({
            "HP:ROOT":   [],
            "HP:ANCHOR": ["HP:ROOT"],
            "HP:LEAF":   ["HP:ANCHOR"],
        })
        assert lift_hp_to_anchors("HP:LEAF", ont, {"HP:ANCHOR"}) == {"HP:ANCHOR": 1}

    def test_leaf_with_multiple_anchor_ancestors(self):
        result = lift_hp_to_anchors("HP:LEAF", ONT, ANCHORS)
        assert set(result.keys()) == ANCHORS

    def test_distances_are_nonneg_int(self):
        result = lift_hp_to_anchors("HP:LEAF", ONT, ANCHORS)
        assert all(isinstance(d, int) and d >= 0 for d in result.values())

    def test_leaf_distances_to_direct_anchors_are_one(self):
        result = lift_hp_to_anchors("HP:LEAF", ONT, ANCHORS)
        # HP:LEAF's direct parents are the two anchors — distance must be 1.
        assert result == {"HP:ANCHOR_A": 1, "HP:ANCHOR_B": 1}

    def test_no_anchor_ancestor_returns_empty(self):
        assert lift_hp_to_anchors("HP:ORPHAN", ONT, ANCHORS) == {}

    def test_empty_anchor_set_returns_empty(self):
        assert lift_hp_to_anchors("HP:LEAF", ONT, set()) == {}


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
        assert list(result.columns) == ["hp_id", "mp_top_id", "anchor_hp_id", "hp_to_anchor_distance"]

    def test_anchor_term_maps_directly_with_distance_zero(self):
        hp2mp = {"HP:ANCHOR_A": ["MP:001"]}
        result = build_hp_to_mp_top(["HP:ANCHOR_A"], ONT, hp2mp)
        assert len(result) == 1
        row = result.iloc[0]
        assert row["hp_id"] == "HP:ANCHOR_A"
        assert row["mp_top_id"] == "MP:001"
        assert row["anchor_hp_id"] == "HP:ANCHOR_A"
        assert row["hp_to_anchor_distance"] == 0

    def test_duplicate_inputs_deduplicated(self):
        hp2mp = {"HP:ANCHOR_A": ["MP:001"]}
        ont = _MockOntology({
            "HP:ROOT":     [],
            "HP:ANCHOR_A": ["HP:ROOT"],
            "HP:LEAF":     ["HP:ANCHOR_A"],
        })
        result = build_hp_to_mp_top(["HP:LEAF", "HP:LEAF"], ont, hp2mp)
        assert len(result) == 1

    def test_output_columns(self):
        result = build_hp_to_mp_top(["HP:LEAF"], ONT, {"HP:ANCHOR_A": ["MP:001"]})
        assert list(result.columns) == ["hp_id", "mp_top_id", "anchor_hp_id", "hp_to_anchor_distance"]

    def test_many_to_many_one_anchor_multiple_mp(self):
        # One HP anchor → multiple MP top terms.
        hp2mp = {"HP:ANCHOR_A": ["MP:001", "MP:002"]}
        result = build_hp_to_mp_top(["HP:ANCHOR_A"], ONT, hp2mp)
        assert set(result["mp_top_id"]) == {"MP:001", "MP:002"}
        assert all(result["anchor_hp_id"] == "HP:ANCHOR_A")
        assert all(result["hp_to_anchor_distance"] == 0)

    def test_min_distance_anchor_is_chosen_when_multiple_paths(self):
        # HP:LEAF → ANCHOR_A (dist 1) and HP:LEAF → ANCHOR_B (dist 1).
        # Both anchors map to the same MP term — the closer anchor wins.
        # Here both are equidistant; verify one is chosen and distance is 1.
        hp2mp = {"HP:ANCHOR_A": ["MP:001"], "HP:ANCHOR_B": ["MP:001"]}
        result = build_hp_to_mp_top(["HP:LEAF"], ONT, hp2mp)
        assert len(result) == 1
        assert result.iloc[0]["mp_top_id"] == "MP:001"
        assert result.iloc[0]["hp_to_anchor_distance"] == 1
