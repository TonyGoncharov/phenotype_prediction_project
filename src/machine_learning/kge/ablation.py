"""
src/machine_learning/kge/ablation.py - layer ablation study for the KGE track.

Trains one KGE model per ablation condition, evaluates AUPRC + N-way G→P
ranking metrics, and writes a comparison table.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd
import torch
from pykeen.triples import TriplesFactory
from src.machine_learning.common.results import (
    comparison_dataframe,
    patch_json,
    save_comparison_csv,
)
from src.machine_learning.kge import conditions as kge_conditions
from src.machine_learning.kge.config import (
    DEFAULT_CONFIG,
    target_relation_for_species,
)
from src.machine_learning.kge.evaluate import evaluate_auprc
from src.machine_learning.kge.loader import build_triples_factory, resolve_biocypher_dir
from src.machine_learning.kge.predict import GenePhenotypePredictor
from src.machine_learning.kge.train import run_training, split_triples_factory

logger = logging.getLogger(__name__)

CONDITIONS = kge_conditions.CONDITIONS
HUMAN_CONDITIONS = kge_conditions.HUMAN_CONDITIONS
MOUSE_CONDITIONS = kge_conditions.MOUSE_CONDITIONS
SPECIES_CONDITIONS = kge_conditions.SPECIES_CONDITIONS

_FAST_PRESET = {
    'embedding_dim': 64,
    'num_epochs': 30,
    'num_negs_per_pos': 32,
    'batch_size': 1024,
}


# ── Contamination-safe split ──────────────────────────────────────────────────


def make_condition_splits(
    data_dir: str | Path,
    include_relations: set[str] | None,
    locked_test_pairs: dict[str, set[str]],
    species: str = 'human',
    target_relation: str | None = None,
    train_frac: float = DEFAULT_CONFIG['train_frac'],
    val_frac: float = DEFAULT_CONFIG['val_frac'],
    seed: int = DEFAULT_CONFIG['random_seed'],
) -> tuple[TriplesFactory, TriplesFactory, TriplesFactory]:
    """
    Return (train_tf, val_tf, test_tf) with the test fold locked to fixed pairs.

    Gene→MP pairs in *locked_test_pairs* are pre-assigned to test and never
    appear in training, regardless of entity-to-integer mapping differences
    across conditions.

    Args:
        data_dir: BioCypher output directory.
        include_relations: Relation set for this condition (None = all).
        locked_test_pairs: mp_term_id → {gene_id, …} from the full-graph split.
        species: 'human' or 'mouse' (selects the correct phenotype relation).
        target_relation: Explicit gene→phenotype target relation.
        train_frac: Training fraction of remaining (non-test) triples.
        val_frac: Validation fraction of remaining triples.
        seed: RNG seed.

    Returns:
        (train_tf, val_tf, test_tf) - all sharing the same entity/relation vocabulary.
    """
    pheno_rel = target_relation_for_species(species, target_relation)
    tf = build_triples_factory(
        data_dir,
        include_relations=include_relations,
        target_relation=pheno_rel,
    )
    rel_id = tf.relation_to_id.get(pheno_rel)

    if rel_id is None:
        logger.warning(
            "'%s' absent from condition - falling back to standard split", pheno_rel
        )
        return split_triples_factory(tf, train_frac, val_frac, seed)

    id_to_entity: dict[int, str] = {v: k for k, v in tf.entity_to_id.items()}

    locked_mask = torch.zeros(tf.num_triples, dtype=torch.bool)
    gp_indices = (tf.mapped_triples[:, 1] == rel_id).nonzero(as_tuple=True)[0]
    for idx in gp_indices.tolist():
        h, _, t = tf.mapped_triples[idx].tolist()
        gene = id_to_entity.get(h, '')
        mp_trm = id_to_entity.get(t, '')
        if mp_trm in locked_test_pairs and gene in locked_test_pairs[mp_trm]:
            locked_mask[idx] = True

    n_locked = int(locked_mask.sum().item())
    logger.info(
        'Condition split: %d / %d G→P triples locked into test fold',
        n_locked,
        int((tf.mapped_triples[:, 1] == rel_id).sum().item()),
    )

    if n_locked == 0:
        logger.warning(
            'No locked test pairs found in this condition - standard split used.'
        )
        return split_triples_factory(tf, train_frac, val_frac, seed)

    locked_triples = tf.mapped_triples[locked_mask]
    remainder_triples = tf.mapped_triples[~locked_mask]

    remainder_tf = TriplesFactory(
        mapped_triples=remainder_triples,
        entity_to_id=tf.entity_to_id,
        relation_to_id=tf.relation_to_id,
    )

    adj_train_ratio = train_frac / (train_frac + val_frac)
    train_tf, val_tf = remainder_tf.split(ratios=[adj_train_ratio], random_state=seed)

    test_tf = TriplesFactory(
        mapped_triples=locked_triples,
        entity_to_id=tf.entity_to_id,
        relation_to_id=tf.relation_to_id,
    )

    logger.info(
        'Contamination-safe split - train: %d  val: %d  test: %d',
        train_tf.num_triples,
        val_tf.num_triples,
        test_tf.num_triples,
    )
    return train_tf, val_tf, test_tf


def get_fixed_test_positives(
    data_dir: str | Path,
    species: str = 'human',
    target_relation: str | None = None,
    train_frac: float = DEFAULT_CONFIG['train_frac'],
    val_frac: float = DEFAULT_CONFIG['val_frac'],
    seed: int = DEFAULT_CONFIG['random_seed'],
) -> dict[str, set[str]]:
    """
    Build the canonical gene→phenotype test set from the full graph.

    Loads all relation types, splits with the given seed, and returns the
    gene→phenotype pairs that landed in the test fold as:
        { mp_term_id → set of gene_id strings }

    Computed once before the ablation loop so every condition evaluates against
    exactly the same positive pairs.

    Args:
        data_dir: BioCypher output directory.
        species: 'human' or 'mouse'.
        target_relation: Explicit gene→phenotype target relation.
        train_frac: Training fraction.
        val_frac: Validation fraction.
        seed: RNG seed.

    Returns:
        Dict mapping MP term ID → set of positive gene IDs in the test fold.
    """
    pheno_rel = target_relation_for_species(species, target_relation)
    logger.info('Computing fixed G→P test positives from full graph …')
    tf = build_triples_factory(
        data_dir,
        include_relations=None,
        target_relation=pheno_rel,
    )
    _, _, test_tf = split_triples_factory(
        tf, train_frac=train_frac, val_frac=val_frac, seed=seed
    )

    rel_id = test_tf.relation_to_id.get(pheno_rel)
    if rel_id is None:
        raise ValueError(
            f"'{pheno_rel}' not found in the full graph. "
            'Cannot build fixed test positives.'
        )

    mask = test_tf.mapped_triples[:, 1] == rel_id
    id_to_entity = {v: k for k, v in test_tf.entity_to_id.items()}

    positives: dict[str, set[str]] = {}
    for head_id, _, tail_id in test_tf.mapped_triples[mask].tolist():
        gene = id_to_entity[head_id]
        mp_term = id_to_entity[tail_id]
        positives.setdefault(mp_term, set()).add(gene)

    total = sum(len(v) for v in positives.values())
    logger.info(
        'Fixed test set: %d gene→MP pairs across %d MP terms',
        total,
        len(positives),
    )
    return positives


# ── Core ablation loop ────────────────────────────────────────────────────────


def run_ablation(
    data_dir: str | Path | None = None,
    out_dir: str | Path | None = None,
    conditions: dict[str, set[str] | None] | None = None,
    species: str = 'human',
    target_relation: str | None = None,
    fast: bool = False,
    target_steps: int | None = None,
    model: str = 'RotatE',
    config: dict | None = None,
) -> pd.DataFrame:
    """
    Train and evaluate one KGE model per ablation condition.

    All conditions share the same fixed gene→phenotype test set (derived from
    the full-graph split) so AUPRC values are directly comparable.

    Note on MRR/Hits@K: these rank against all entities in each condition's
    graph, so the candidate pool differs across conditions (removing PPI
    shrinks entities from ~60k to ~36k, artificially boosting ranking metrics).
    MRR/Hits@K are reported for completeness but AUPRC is the primary metric.

    Args:
        data_dir: BioCypher output directory. None → auto-detect.
        out_dir: Root output directory. Condition results go in subdirs.
        conditions: Dict of condition_name → relation set. None → species default.
        species: 'human' or 'mouse'.
        target_relation: Explicit gene→phenotype target relation.
        fast: Smoke-test preset (dim=64, 30 epochs).
        target_steps: Gradient-step budget per condition.
        model: KGE model name.
        config: Optional training/evaluation settings merged over DEFAULT_CONFIG.

    Returns:
        DataFrame with one row per condition, written to <out_dir>/comparison_summary.csv.
    """
    base_cfg = {**DEFAULT_CONFIG, **(config or {})}
    species = species.lower()
    target_relation = target_relation_for_species(
        species,
        target_relation or base_cfg.get('target_relation'),
    )
    data_dir = resolve_biocypher_dir(data_dir, species=species)
    conditions = conditions or SPECIES_CONDITIONS.get(species, HUMAN_CONDITIONS)
    if out_dir is None:
        out_dir = f'results/ml/{species}/kge/{model.lower()}'
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    fixed_positives = get_fixed_test_positives(
        data_dir,
        species=species,
        target_relation=target_relation,
        train_frac=base_cfg['train_frac'],
        val_frac=base_cfg['val_frac'],
        seed=base_cfg['random_seed'],
    )

    rows: list[dict] = []

    for name, rels in conditions.items():
        cond_dir = out_path / name
        rels_desc = 'all' if rels is None else '|'.join(sorted(rels))

        logger.info('=' * 70)
        logger.info('Condition: %s', name)
        logger.info('Relations: %s', rels_desc)
        logger.info('=' * 70)

        try:
            cfg = {
                **base_cfg,
                # Disable early stopping: ES triggers full-entity evaluation every
                # es_frequency epochs. Budget is controlled by target_steps, so ES
                # adds no benefit and doubles peak RAM usage here.
                'es_patience': base_cfg.get('es_patience', 0),
                **(fast and _FAST_PRESET or {}),
                'model': model,
                'include_relations': rels,
                **(({'target_steps': target_steps}) if target_steps else {}),
            }

            train_tf, val_tf, test_tf = make_condition_splits(
                data_dir=data_dir,
                include_relations=rels,
                locked_test_pairs=fixed_positives,
                species=species,
                target_relation=target_relation,
                train_frac=base_cfg['train_frac'],
                val_frac=base_cfg['val_frac'],
                seed=base_cfg['random_seed'],
            )

            summary = run_training(
                data_dir=str(data_dir),
                out_dir=cond_dir,
                config=cfg,
                species=species,
                target_relation=target_relation,
                train_tf=train_tf,
                val_tf=val_tf,
                test_tf=test_tf,
            )

            auprc_result = _run_auprc(
                cond_dir,
                fixed_positives=fixed_positives,
                target_relation=target_relation,
            )

        except Exception:
            logger.exception("Condition '%s' failed - skipping", name)
            rows.append({'condition': name, 'error': 'failed'})
            continue

        m = summary.get('gp_metrics') or summary.get('metrics', {})
        rows.append(
            {
                'condition': name,
                'species': species,
                'target_relation': target_relation,
                'relations': rels_desc,
                'mean_auprc': auprc_result.get('mean_auprc'),
                'auprc_classes': auprc_result.get('num_classes'),
                'mrr': m.get('gp_mrr', m.get('mrr')),
                'hits_at_1': m.get('gp_hits_at_1', m.get('hits_at_1')),
                'hits_at_3': m.get('gp_hits_at_3', m.get('hits_at_3')),
                'hits_at_10': m.get('gp_hits_at_10', m.get('hits_at_10')),
                'gp_n_candidates': m.get('gp_n_candidates'),
                'gp_n_test': m.get('gp_n_test'),
            }
        )
        logger.info(
            'Done - AUPRC=%.4f  MRR=%.4f  Hits@10=%.4f',
            auprc_result.get('mean_auprc') or 0,
            m.get('gp_mrr', m.get('mrr')) or 0,
            m.get('gp_hits_at_10', m.get('hits_at_10')) or 0,
        )

    df = comparison_dataframe(rows)
    save_comparison_csv(df, out_path)
    logger.info('Comparison summary → %s', out_path / 'comparison_summary.csv')
    return df


# ── Helpers ───────────────────────────────────────────────────────────────────


def _run_auprc(
    cond_dir: Path,
    fixed_positives: dict[str, set[str]] | None = None,
    target_relation: str | None = None,
) -> dict:
    """
    Load test_triples, compute AUPRC, write per-class CSV, patch summary.json.
    """
    test_tf_path = cond_dir / 'test_triples'
    if not test_tf_path.exists():
        logger.warning('test_triples/ not found in %s - AUPRC skipped', cond_dir)
        return {'mean_auprc': None, 'num_classes': 0}

    predictor = GenePhenotypePredictor.from_directory(
        cond_dir,
        relation=target_relation,
    )
    test_tf = TriplesFactory.from_path_binary(test_tf_path)
    result = evaluate_auprc(
        predictor,
        test_tf,
        relation=predictor.relation,
        override_positives=fixed_positives,
    )

    result['per_class'].to_csv(cond_dir / 'auprc_per_phenotype.csv', index=False)

    sp = cond_dir / 'summary.json'
    if sp.exists():
        patch_json(
            sp,
            {
                'metrics': {
                    **json.loads(sp.read_text()).get('metrics', {}),
                    'mean_auprc': result['mean_auprc'],
                    'auprc_classes': result['num_classes'],
                }
            },
        )

    return result
