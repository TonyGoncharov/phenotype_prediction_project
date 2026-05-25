"""
src/machine_learning/kge/conditions.py - dependency-free KGE ablation condition definitions.

This module intentionally avoids importing PyKEEN so common tests and lightweight
utilities can inspect condition names without requiring the KGE training stack.
"""

from __future__ import annotations

from src.machine_learning.kge.config import (
    GENE_PHENOTYPE_RELATION_HUMAN,
    GENE_PHENOTYPE_RELATION_MOUSE,
)

# Relation groups -- human
_H_PHENO = {GENE_PHENOTYPE_RELATION_HUMAN}
_H_GO = {'HumanGeneHasGoTerm'}
_H_EXPR = {'HumanGeneExpressedInTissue', 'TissueMappedToUberon'}
_H_PPI = {'HumanGeneEncodesProtein', 'HumanProteinInteractsWith'}

# Relation groups -- mouse
_M_PHENO = {GENE_PHENOTYPE_RELATION_MOUSE}
_M_GO = {'MouseGeneHasGoTerm'}
_M_EXPR = {'MouseGeneExpressedInTissue'}
_M_PPI = {'MouseGeneEncodesProtein', 'MouseProteinInteractsWith'}

HUMAN_CONDITIONS: dict[str, set[str] | None] = {
    'all': None,
    'no_go': _H_PHENO | _H_EXPR | _H_PPI,
    'no_ppi': _H_PHENO | _H_GO | _H_EXPR,
    'no_expr': _H_PHENO | _H_GO | _H_PPI,
    'pheno_only': _H_PHENO,
}

MOUSE_CONDITIONS: dict[str, set[str] | None] = {
    'all': None,
    'no_go': _M_PHENO | _M_EXPR | _M_PPI,
    'no_ppi': _M_PHENO | _M_GO | _M_EXPR,
    'no_expr': _M_PHENO | _M_GO | _M_PPI,
    'pheno_only': _M_PHENO,
}

# Backward-compatible default.
CONDITIONS = HUMAN_CONDITIONS

SPECIES_CONDITIONS = {'human': HUMAN_CONDITIONS, 'mouse': MOUSE_CONDITIONS}
