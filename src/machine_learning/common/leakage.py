"""
src/machine_learning/common/leakage.py - leakage guards for ML experiments.

These checks intentionally fail early on direct target leakage and emit warnings
for suspicious relation states that deserve attention but are not, by
themselves, proof that target labels reached the model.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import pandas as pd
from src.machine_learning.kge.config import (
    SPECIES_TARGET_RELATIONS,
    target_relation_for_species,
)

if TYPE_CHECKING:
    from pykeen.triples import TriplesFactory

logger = logging.getLogger(__name__)


def check_kge_target_leakage(
    train_tf: TriplesFactory,
    test_tf: TriplesFactory,
    species: str = 'human',
    target_relation: str | None = None,
    val_tf: TriplesFactory | None = None,
) -> None:
    """
    Ensure test target gene→phenotype triples are absent from KGE training triples.

    Args:
        train_tf: Training triples.
        test_tf: Test triples.
        species: Species used to resolve the target relation.
        target_relation: Explicit target relation; overrides species default.
        val_tf: Optional validation triples, checked only for suspicious overlap.

    Raises:
        ValueError: If any exact target triple from test appears in training.
    """
    target_relation = target_relation_for_species(species, target_relation)
    train_targets = _target_pairs_from_tf(
        train_tf, target_relation, split_name='training'
    )
    test_targets = _target_pairs_from_tf(test_tf, target_relation, split_name='test')

    if not test_targets:
        logger.warning(
            "KGE leakage check: no test triples found for target relation '%s' "
            '(species=%s). Gene→phenotype evaluation may be empty.',
            target_relation,
            species,
        )
        return

    if not train_targets:
        logger.warning(
            "KGE leakage check: target relation '%s' is absent from training triples "
            '(species=%s). This may be intentional for an ablation, but candidate '
            'phenotype coverage can be weak.',
            target_relation,
            species,
        )

    leaked = train_targets & test_targets
    if leaked:
        examples = _format_pair_examples(leaked)
        raise ValueError(
            'Direct KGE target leakage detected: '
            f'{len(leaked)} test {target_relation} gene→phenotype triple(s) are '
            f'also present in training triples for species={species}. '
            f'Examples: {examples}'
        )

    if val_tf is not None:
        val_overlap = (
            _target_pairs_from_tf(
                val_tf,
                target_relation,
                split_name='validation',
                warn_if_absent=False,
            )
            & test_targets
        )
        if val_overlap:
            logger.warning(
                'KGE leakage check: %d test target triple(s) also appear in validation. '
                'They are not in training, but early stopping/evaluation may be optimistic. '
                'Examples: %s',
                len(val_overlap),
                _format_pair_examples(val_overlap),
            )

    reversed_overlap = {(tail, head) for head, tail in test_targets} & train_targets
    if reversed_overlap:
        logger.warning(
            'KGE leakage check: %d reversed target pair(s) from test appear in training '
            "under relation '%s'. This is not exact leakage but is suspicious. Examples: %s",
            len(reversed_overlap),
            target_relation,
            _format_pair_examples(reversed_overlap),
        )

    logger.info(
        "KGE leakage check passed for '%s' (species=%s): %d train target pairs, "
        '%d test target pairs, 0 direct overlaps.',
        target_relation,
        species,
        len(train_targets),
        len(test_targets),
    )


def check_rgcn_context_leakage(
    context_df: pd.DataFrame,
    species: str = 'human',
    target_relation: str | None = None,
    active_rels: set[str] | None = None,
) -> None:
    """
    Ensure R-GCN message-passing context does not contain target phenotype edges.

    Args:
        context_df: Final context graph edges used for message passing.
        species: Species used to resolve the target relation.
        target_relation: Explicit target relation; overrides species default.
        active_rels: Optional configured relation set for the condition.

    Raises:
        ValueError: If the final context graph contains the species target relation.
    """
    target_relation = target_relation_for_species(species, target_relation)
    active_rels = set(active_rels or set())

    if target_relation in active_rels:
        logger.warning(
            "R-GCN leakage check: target relation '%s' was requested in active_rels "
            '(species=%s). The graph builder should still exclude it from message '
            'passing; verify the condition definition.',
            target_relation,
            species,
        )

    other_target_rels = {
        rel
        for sp, rel in SPECIES_TARGET_RELATIONS.items()
        if sp != species and rel in active_rels
    }
    if other_target_rels:
        logger.warning(
            'R-GCN leakage check: active_rels contains phenotype relation(s) for a '
            'different species: %s. This is suspicious for species=%s.',
            sorted(other_target_rels),
            species,
        )

    if context_df.empty:
        logger.warning(
            'R-GCN leakage check: context graph is empty for species=%s. No direct '
            'target leakage is present, but training will have no message-passing edges.',
            species,
        )
        return

    direct = context_df[context_df['relation'] == target_relation]
    if not direct.empty:
        examples = _format_df_examples(direct)
        raise ValueError(
            'Direct R-GCN target leakage detected: '
            f'{len(direct)} {target_relation} gene→phenotype edge(s) are present '
            f'in the message-passing context graph for species={species}. '
            f'Phenotype labels must be used only downstream. Examples: {examples}'
        )

    other_context_targets = sorted(
        rel
        for rel in SPECIES_TARGET_RELATIONS.values()
        if rel != target_relation and rel in set(context_df['relation'].unique())
    )
    if other_context_targets:
        logger.warning(
            'R-GCN leakage check: context graph contains phenotype-like target '
            'relation(s) for another species: %s. This is not direct leakage for '
            'species=%s, but it is suspicious.',
            other_context_targets,
            species,
        )

    phenotype_named = sorted(
        rel
        for rel in context_df['relation'].unique()
        if rel != target_relation and ('Phenotype' in rel or 'HasMpTopTerm' in rel)
    )
    if phenotype_named:
        logger.warning(
            'R-GCN leakage check: context graph contains phenotype-named relation(s) '
            'that are not the configured target relation %r: %s.',
            target_relation,
            phenotype_named,
        )

    logger.info(
        "R-GCN leakage check passed for '%s' (species=%s): %d context edges, "
        '0 direct target edges.',
        target_relation,
        species,
        len(context_df),
    )


def _target_pairs_from_tf(
    tf: TriplesFactory,
    target_relation: str,
    split_name: str,
    warn_if_absent: bool = True,
) -> set[tuple[str, str]]:
    """
    Extract (gene_id, mp_term_id) pairs for a given relation from a TriplesFactory.
    """
    rel_id = tf.relation_to_id.get(target_relation)
    if rel_id is None:
        if warn_if_absent:
            logger.warning(
                "KGE leakage check: target relation '%s' absent from %s triples.",
                target_relation,
                split_name,
            )
        return set()

    mask = tf.mapped_triples[:, 1] == rel_id
    triples = tf.mapped_triples[mask]
    id_to_entity = {idx: entity for entity, idx in tf.entity_to_id.items()}
    return {
        (id_to_entity[int(head_id)], id_to_entity[int(tail_id)])
        for head_id, _, tail_id in triples.tolist()
    }


def _format_pair_examples(pairs: set[tuple[str, str]], limit: int = 5) -> str:
    """
    Format up to *limit* (head, tail) pairs as a semicolon-separated string.
    """
    examples = sorted(pairs)[:limit]
    return '; '.join(f'{head} -> {tail}' for head, tail in examples)


def _format_df_examples(df: pd.DataFrame, limit: int = 5) -> str:
    """
    Format up to *limit* rows from a triples DataFrame as a semicolon-separated string.
    """
    examples = df[['head', 'relation', 'tail']].head(limit).itertuples(index=False)
    return '; '.join(
        f'{head} -[{relation}]-> {tail}' for head, relation, tail in examples
    )
