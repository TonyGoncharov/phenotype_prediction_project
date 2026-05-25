"""
src/machine_learning/rgcn/graph_builder.py - heterogeneous context graph construction for R-GCN.

Builds the PyG-compatible edge tensors that the R-GCN encoder trains on.

Critical invariant
------------------
Phenotype edges (HumanGeneHasMpTopTerm / MouseGeneHasMpTopTerm) are NEVER
included in the context graph. They are separated into pheno_df for use as
classification labels only. Including them here would leak target associations
into the encoder and inflate downstream AUPRC.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import pandas as pd
import torch
from src.machine_learning.common.leakage import check_rgcn_context_leakage
from src.machine_learning.rgcn.config import (
    PHENOTYPE_REL_HUMAN,  # noqa: F401 - re-exported for test and external callers
    PHENOTYPE_REL_MOUSE,  # noqa: F401 - re-exported for test and external callers
    target_relation_for_species,
)

logger = logging.getLogger(__name__)

# ── Relation name constants ───────────────────────────────────────────────────

# Full context relation sets per species (phenotype excluded).
CONTEXT_RELS_HUMAN: set[str] = {
    'HumanGeneEncodesProtein',
    'HumanGeneExpressedInTissue',
    'HumanGeneHasGoTerm',
    'HumanProteinInteractsWith',
    'TissueMappedToUberon',
}

CONTEXT_RELS_MOUSE: set[str] = {
    'MouseGeneEncodesProtein',
    'MouseGeneExpressedInTissue',
    'MouseGeneHasGoTerm',
    'MouseProteinInteractsWith',
}

_SPECIES_CONTEXT_RELS = {'human': CONTEXT_RELS_HUMAN, 'mouse': CONTEXT_RELS_MOUSE}
_SPECIES_ENCODES_REL = {
    'human': 'HumanGeneEncodesProtein',
    'mouse': 'MouseGeneEncodesProtein',
}
_SPECIES_PPI_REL = {
    'human': 'HumanProteinInteractsWith',
    'mouse': 'MouseProteinInteractsWith',
}
_SPECIES_EXPR_REL = {
    'human': 'HumanGeneExpressedInTissue',
    'mouse': 'MouseGeneExpressedInTissue',
}
_SPECIES_GO_REL = {
    'human': 'HumanGeneHasGoTerm',
    'mouse': 'MouseGeneHasGoTerm',
}
_SPECIES_UBERON_REL = {
    'human': 'TissueMappedToUberon',
    'mouse': None,
}


@dataclass
class ContextGraph:
    """
    Encapsulates the tensors and mappings needed to train an R-GCN encoder.

    Attributes:
        edge_index: Bidirectional edge connectivity, shape (2, 2 * n_edges).
        edge_type: Relation type indices (forward + inverse), shape (2 * n_edges,).
        num_entities: Total entity count in the context graph.
        num_relations: Total relation types including inverses.
        ent2id: Entity string → integer index.
        rel2id: Relation string → integer index (forward relations only).
        context_df: DataFrame of context edges (head, relation, tail).
        pheno_df: DataFrame of phenotype edges used as labels (head, relation, tail).
        all_genes: Set of gene IDs that have at least one phenotype annotation.
    """

    edge_index: torch.Tensor
    edge_type: torch.Tensor
    num_entities: int
    num_relations: int
    ent2id: dict[str, int]
    rel2id: dict[str, int]
    context_df: pd.DataFrame
    pheno_df: pd.DataFrame
    all_genes: set[str]


def build_context_graph(
    all_df: pd.DataFrame,
    species: str = 'human',
    active_rels: set[str] | None = None,
    target_relation: str | None = None,
) -> ContextGraph:
    """
    Build the R-GCN context graph from a full BioCypher triples DataFrame.

    Phenotype edges are separated from context edges and returned in pheno_df.
    The context graph is always bidirectional (forward + inverse relation types).

    Args:
        all_df: DataFrame with columns [head, relation, tail] containing all
            triples loaded from BioCypher (including phenotype edges).
        species: 'human' or 'mouse'. Determines phenotype and context relation names.
        active_rels: Subset of context relations to include. None uses the full
            species-default context set (CONTEXT_RELS_HUMAN or CONTEXT_RELS_MOUSE).
        target_relation: Optional phenotype relation override.

    Returns:
        ContextGraph dataclass with edge tensors, entity/relation mappings,
        and the phenotype DataFrame.
    """
    pheno_rel = target_relation_for_species(species, target_relation)
    default_rels = _SPECIES_CONTEXT_RELS.get(species, CONTEXT_RELS_HUMAN)
    used_rels = active_rels if active_rels is not None else default_rels

    encodes_rel = _SPECIES_ENCODES_REL.get(species)
    ppi_rel = _SPECIES_PPI_REL.get(species)
    expr_rel = _SPECIES_EXPR_REL.get(species)
    go_rel = _SPECIES_GO_REL.get(species)
    uberon_rel = _SPECIES_UBERON_REL.get(species)

    # Phenotype edges → labels only, never in the graph.
    pheno_df = all_df[all_df['relation'] == pheno_rel].copy()
    all_genes = set(pheno_df['head'].unique())

    logger.info('Genes with phenotype labels: %d', len(all_genes))
    logger.info('Gene→Phenotype triples: %d', len(pheno_df))
    logger.info('Unique MP top-terms: %d', pheno_df['tail'].nunique())

    # ── Build context sub-graphs ──────────────────────────────────────────────
    ctx = all_df[all_df['relation'].isin(used_rels)]

    enc = _filter_edges(ctx, encodes_rel)
    enc = enc[enc['head'].isin(all_genes)] if len(enc) else enc
    proteins = set(enc['tail'].unique())

    ppi = _filter_edges(ctx, ppi_rel)
    ppi = (
        ppi[ppi['head'].isin(proteins) | ppi['tail'].isin(proteins)]
        if len(ppi)
        else ppi
    )

    expr = _filter_edges(ctx, expr_rel)
    expr = expr[expr['head'].isin(all_genes)] if len(expr) else expr
    tissues = set(expr['tail'].unique())

    go = _filter_edges(ctx, go_rel)
    go = go[go['head'].isin(all_genes)] if len(go) else go

    uberon_parts = []
    if uberon_rel:
        ub = _filter_edges(ctx, uberon_rel)
        ub = ub[ub['head'].isin(tissues)] if len(ub) else ub
        uberon_parts.append(ub)

    context_df = pd.concat(
        [enc, ppi, expr, go] + uberon_parts, ignore_index=True
    ).drop_duplicates()
    check_rgcn_context_leakage(
        context_df=context_df,
        species=species,
        target_relation=pheno_rel,
        active_rels=used_rels,
    )

    logger.info('R-GCN context graph (leakage-free - no phenotype edges):')
    logger.info(
        '  Entities: %d', len(set(context_df['head']) | set(context_df['tail']))
    )
    logger.info(
        '  Encodes: %d | PPI: %d | Expression: %d | GO: %d | Uberon: %d',
        len(enc),
        len(ppi),
        len(expr),
        len(go),
        len(uberon_parts[0]) if uberon_parts else 0,
    )

    # ── Build PyG tensors ─────────────────────────────────────────────────────
    entities = sorted(set(context_df['head']) | set(context_df['tail']))
    relations = sorted(context_df['relation'].unique())
    ent2id = {e: i for i, e in enumerate(entities)}
    rel2id = {r: i for i, r in enumerate(relations)}

    heads = torch.tensor([ent2id[h] for h in context_df['head']], dtype=torch.long)
    tails = torch.tensor([ent2id[t] for t in context_df['tail']], dtype=torch.long)
    rtypes = torch.tensor([rel2id[r] for r in context_df['relation']], dtype=torch.long)

    # Bidirectional: forward + inverse (standard R-GCN practice).
    edge_index = torch.stack(
        [
            torch.cat([heads, tails]),
            torch.cat([tails, heads]),
        ]
    )
    edge_type = torch.cat([rtypes, rtypes + len(relations)])
    num_relations_total = len(relations) * 2

    n_missing = sum(1 for g in all_genes if g not in ent2id)
    if n_missing:
        logger.warning(
            '%d genes have phenotype labels but no context edges - '
            'they will receive zero embeddings (not dropped).',
            n_missing,
        )

    logger.info(
        'Context graph ready - %d entities | %d relation types (×2 = %d with inverses) | %d edges',
        len(entities),
        len(relations),
        num_relations_total,
        len(context_df),
    )

    return ContextGraph(
        edge_index=edge_index,
        edge_type=edge_type,
        num_entities=len(entities),
        num_relations=num_relations_total,
        ent2id=ent2id,
        rel2id=rel2id,
        context_df=context_df,
        pheno_df=pheno_df,
        all_genes=all_genes,
    )


def _filter_edges(df: pd.DataFrame, relation: str | None) -> pd.DataFrame:
    """
    Return rows matching *relation*, or an empty DataFrame if relation is None.
    """
    if relation is None or relation not in df['relation'].values:
        return pd.DataFrame(columns=['head', 'relation', 'tail'])
    return df[df['relation'] == relation].copy()
