"""Tests for src/layers/ppi_export.py

Covers _detect_column_map() — the BioGRID column-format auto-detection
logic that distinguishes v4 (pre-5.0) and v5 (5.0+) file schemas.
"""

import pytest

from src.layers.ppi_export import _BIOGRID_FIELD_ALIASES, _detect_column_map

# Build canonical column lists directly from the aliases table so these
# tests stay in sync if new fields are ever added to _BIOGRID_FIELD_ALIASES.
_V4_COLUMNS = [aliases[0] for aliases in _BIOGRID_FIELD_ALIASES.values()]
_V5_COLUMNS = [aliases[1] for aliases in _BIOGRID_FIELD_ALIASES.values()]


class TestDetectColumnMap:
    @pytest.mark.parametrize("columns,expected_raw", [
        (_V4_COLUMNS, "Official Symbol Interactor A"),
        (_V5_COLUMNS, "Gene Name Interactor A"),
    ])
    def test_known_format_returns_full_mapping(self, columns, expected_raw):
        result = _detect_column_map(columns)
        assert len(result) == len(_BIOGRID_FIELD_ALIASES)
        assert result[expected_raw] == "symbol_a"

    def test_v4_maps_symbol_columns(self):
        result = _detect_column_map(_V4_COLUMNS)
        assert result["Official Symbol Interactor A"] == "symbol_a"
        assert result["Official Symbol Interactor B"] == "symbol_b"

    def test_v5_maps_symbol_columns(self):
        result = _detect_column_map(_V5_COLUMNS)
        assert result["Gene Name Interactor A"] == "symbol_a"
        assert result["Gene Name Interactor B"] == "symbol_b"

    def test_extra_columns_are_ignored(self):
        cols = _V4_COLUMNS + ["ExtraCol1", "ExtraCol2"]
        result = _detect_column_map(cols)
        assert len(result) == len(_BIOGRID_FIELD_ALIASES)

    def test_unknown_columns_raise_value_error(self):
        with pytest.raises(ValueError, match="Cannot map"):
            _detect_column_map(["ColA", "ColB", "ColC"])

    def test_missing_one_required_field_raises(self):
        # Drop "Official Symbol Interactor A" — v4 no longer fully matches,
        # v5 also won't match (different symbol column names), so ValueError.
        incomplete = [c for c in _V4_COLUMNS if c != "Official Symbol Interactor A"]
        with pytest.raises(ValueError):
            _detect_column_map(incomplete)

    def test_error_message_names_missing_fields(self):
        with pytest.raises(ValueError, match="symbol_a"):
            incomplete = [c for c in _V4_COLUMNS if c != "Official Symbol Interactor A"]
            _detect_column_map(incomplete)

    def test_empty_column_list_raises(self):
        with pytest.raises(ValueError):
            _detect_column_map([])
