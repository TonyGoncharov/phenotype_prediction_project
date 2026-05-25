"""
src/knowledge_graph/layers/_mp_ontology.py - shared MP ontology helpers.

Shared constants and traversal functions for the Mammalian Phenotype (MP)
ontology used by both human and mouse gene-to-phenotype export pipelines.
"""

from __future__ import annotations

MP_ROOT = 'MP:0000001'


def top_children(ont, root_id: str) -> set[str]:
    """
    Return the direct children of root_id in the MP ontology.

    Args:
        ont: pronto.Ontology instance.
        root_id: Root MP term ID.

    Returns:
        Set of direct child MP term IDs.
    """
    root = ont[root_id]
    tops = {t.id for t in root.subclasses(distance=1)}
    tops.discard(root_id)
    return tops


def collapse_to_top(ont, term_id: str, top_set: set[str], root_id: str) -> list[str]:
    """
    Map term_id to any top-level MP terms it descends from.

    Args:
        ont: pronto.Ontology instance.
        term_id: MP term ID to map.
        top_set: Set of top-level MP term IDs.
        root_id: Root MP term ID to exclude from results.

    Returns:
        Sorted list of top-level MP term IDs, or empty list if absent.
    """
    if term_id not in ont:
        return []
    anc = {t.id for t in ont[term_id].superclasses()}
    anc.add(term_id)
    return sorted((anc & top_set) - {root_id})
