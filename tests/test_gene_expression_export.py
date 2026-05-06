"""Tests for src/layers/gene_expression_export.py

Covers _tissue_id() — the normalisation function that converts GTEx
tissue display names into stable graph node IDs (GTEX: prefix).
"""

import pytest

from src.layers.gene_expression_export import _tissue_id


@pytest.mark.parametrize("name,expected", [
    # Standard dash-separated GTEx names
    ("Brain - Amygdala",                                "GTEX:Brain_Amygdala"),
    ("Adipose - Subcutaneous",                          "GTEX:Adipose_Subcutaneous"),
    ("Heart - Left Ventricle",                          "GTEX:Heart_Left_Ventricle"),
    ("Lung",                                            "GTEX:Lung"),
    # Parentheses (non-alphanumeric) collapsed into underscores then stripped
    ("Skin - Sun Exposed (Lower leg)",                  "GTEX:Skin_Sun_Exposed_Lower_leg"),
    ("Brain - Anterior cingulate cortex (BA24)",        "GTEX:Brain_Anterior_cingulate_cortex_BA24"),
    ("Brain - Caudate (basal ganglia)",                 "GTEX:Brain_Caudate_basal_ganglia"),
    # Multiple consecutive separators collapse to a single underscore
    ("A  -  B",                                         "GTEX:A_B"),
])
def test_tissue_id_known_names(name, expected):
    assert _tissue_id(name) == expected


def test_tissue_id_always_has_gtex_prefix():
    for name in ["Lung", "Liver", "Heart - Left Ventricle"]:
        assert _tissue_id(name).startswith("GTEX:")


def test_tissue_id_no_leading_or_trailing_underscore():
    result = _tissue_id("(Leading and trailing)")
    suffix = result[len("GTEX:"):]
    assert not suffix.startswith("_")
    assert not suffix.endswith("_")


def test_tissue_id_no_spaces_in_output():
    result = _tissue_id("Brain - Amygdala")
    assert " " not in result


def test_tissue_id_idempotent_on_clean_name():
    # A name with only alphanumeric chars should survive unchanged
    assert _tissue_id("Lung") == "GTEX:Lung"
