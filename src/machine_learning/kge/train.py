"""
src/machine_learning/kge/train.py - PyKEEN training pipeline for the KGE track.

Trains a KGE model (RotatE, ComplEx, DistMult, TransE) on a BioCypher knowledge
graph, evaluates filtered ranking metrics, and saves model artifacts.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import torch
from pykeen.evaluation import RankBasedEvaluator
from pykeen.pipeline import pipeline
from pykeen.triples import TriplesFactory
from src.machine_learning.common.leakage import check_kge_target_leakage
from src.machine_learning.kge.config import DEFAULT_CONFIG, target_relation_for_species
from src.machine_learning.kge.loader import build_triples_factory

logger = logging.getLogger(__name__)

# Models that operate in complex embedding space and are known to fail on MPS.
_COMPLEX_SPACE_MODELS = {'ComplEx', 'RotatE', 'QuatE', 'RotatH', 'OTE'}

# Models known to produce silent NaN/zero gradients on MPS.
_MPS_UNSUPPORTED_MODELS = _COMPLEX_SPACE_MODELS | {'DistMult'}

# Per-model loss overrides for sLCWA training.
# DistMult defaults to MarginRankingLoss(margin=1.0) in PyKEEN. With sLCWA the
# random dot-product scores start near 0, so loss=max(0,1-0)=1.0 from epoch 1
# and stays there. SoftplusLoss has no margin threshold and provides smooth
# gradients at any score magnitude.
_MODEL_LOSS: dict[str, str] = {
    'DistMult': 'softplus',
}

# ── Gradient-step budget ──────────────────────────────────────────────────────


def compute_epochs(
    n_train_triples: int,
    batch_size: int,
    target_steps: int,
) -> int:
    """
    Return the number of epochs that produces approximately target_steps gradient updates.

    Different ablation conditions have different numbers of training triples, so
    a fixed num_epochs gives conditions very different optimisation budgets. This
    function normalises the budget so every condition trains for the same number
    of parameter updates regardless of graph size.

    Args:
        n_train_triples: Number of training triples in this condition.
        batch_size: Mini-batch size.
        target_steps: Desired total gradient updates.

    Returns:
        Number of epochs (at least 1).

    Example:
        >>> compute_epochs(756_000, 2048, 80_000)
        217
    """
    if n_train_triples <= 0 or batch_size <= 0 or target_steps <= 0:
        raise ValueError(
            f'All arguments must be positive integers '
            f'(got n_train_triples={n_train_triples}, '
            f'batch_size={batch_size}, target_steps={target_steps})'
        )
    steps_per_epoch = n_train_triples / batch_size
    return max(1, round(target_steps / steps_per_epoch))


# ── Data splitting ────────────────────────────────────────────────────────────


def split_triples_factory(
    tf: TriplesFactory,
    train_frac: float = DEFAULT_CONFIG['train_frac'],
    val_frac: float = DEFAULT_CONFIG['val_frac'],
    seed: int = DEFAULT_CONFIG['random_seed'],
) -> tuple[TriplesFactory, TriplesFactory, TriplesFactory]:
    """
    Split *tf* into (train, validation, test) TriplesFactory objects.

    PyKEEN's native split guarantees that every entity and relation that
    appears in val/test also appears in training (no cold-start entities).

    Args:
        tf: Source TriplesFactory to split.
        train_frac: Fraction for training.
        val_frac: Fraction for validation.
        seed: RNG seed.

    Returns:
        (train_tf, val_tf, test_tf)
    """
    test_frac = round(1.0 - train_frac - val_frac, 6)
    if test_frac <= 0:
        raise ValueError(
            f'train_frac ({train_frac}) + val_frac ({val_frac}) must be < 1.0'
        )

    train_tf, val_tf, test_tf = tf.split(
        ratios=[train_frac, val_frac, test_frac],
        random_state=seed,
    )

    logger.info(
        'Split - train: %d  |  val: %d  |  test: %d',
        train_tf.num_triples,
        val_tf.num_triples,
        test_tf.num_triples,
    )
    return train_tf, val_tf, test_tf


# ── Training pipeline ─────────────────────────────────────────────────────────


def run_training(
    data_dir: str | Path,
    out_dir: str | Path,
    config: dict | None = None,
    species: str = 'human',
    target_relation: str | None = None,
    device: str = 'auto',
    cache_dir: str | Path | None = None,
    train_tf: TriplesFactory | None = None,
    val_tf: TriplesFactory | None = None,
    test_tf: TriplesFactory | None = None,
) -> dict:
    """
    Full KGE training + evaluation run.

    Returns a summary dict (config + MRR/Hits@K metrics) also written to
    <out_dir>/summary.json.

    Pre-split mode
    --------------
    Pass *train_tf*, *val_tf*, and *test_tf* together to bypass the internal
    load-and-split step. Used by the ablation study to guarantee every condition
    evaluates against the same fixed test positives (no contamination from
    per-condition random splits).

    Args:
        data_dir: BioCypher output directory (used when pre-split is False).
        out_dir: Directory for model artifacts, logs, and summary.
        config: Override dict merged over DEFAULT_CONFIG.
        species: Species label used to resolve the target relation.
        target_relation: Explicit gene→phenotype target relation.
        device: 'auto', 'cpu', 'cuda', or 'mps'.
        cache_dir: TriplesFactory binary cache directory.
        train_tf: Pre-computed training TriplesFactory (all three required together).
        val_tf: Pre-computed validation TriplesFactory.
        test_tf: Pre-computed test TriplesFactory.

    Returns:
        Summary dict written to <out_dir>/summary.json.
    """
    cfg = {**DEFAULT_CONFIG, **(config or {})}
    include_relations = cfg.pop('include_relations', None)
    species = str(cfg.get('species', species)).lower()
    target_relation = target_relation_for_species(
        species=species,
        target_relation=target_relation or cfg.get('target_relation'),
    )
    cfg['species'] = species
    cfg['target_relation'] = target_relation
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pre_split = train_tf is not None and val_tf is not None and test_tf is not None

    # ── 1. Load and split data ────────────────────────────────────────────────
    logger.info('=' * 70)
    if pre_split:
        logger.info('Step 1/3 - Using pre-computed splits (contamination-safe)')
        logger.info(
            '  train: %d  val: %d  test: %d (locked G→P pairs)',
            train_tf.num_triples,
            val_tf.num_triples,
            test_tf.num_triples,
        )
        _num_entities = train_tf.num_entities
        _num_relations = train_tf.num_relations
        _num_triples = train_tf.num_triples + val_tf.num_triples + test_tf.num_triples
    else:
        logger.info('Step 1/3 - Loading triples from %s', data_dir)
        tf = build_triples_factory(
            data_dir,
            include_relations=include_relations,
            cache_dir=cache_dir,
            target_relation=target_relation,
        )
        train_tf, val_tf, test_tf = split_triples_factory(
            tf,
            train_frac=cfg['train_frac'],
            val_frac=cfg['val_frac'],
            seed=cfg['random_seed'],
        )
        _num_entities = tf.num_entities
        _num_relations = tf.num_relations
        _num_triples = tf.num_triples

    _report_relation_split(train_tf, val_tf, test_tf, target_relation=target_relation)
    check_kge_target_leakage(
        train_tf=train_tf,
        val_tf=val_tf,
        test_tf=test_tf,
        species=species,
        target_relation=target_relation,
    )

    # ── 1b. Resolve num_epochs from target_steps ──────────────────────────────
    _explicit_target = (config or {}).get('target_steps')
    if _explicit_target or cfg.get('num_epochs') is None:
        _target = _explicit_target or cfg.get('target_steps')
        if _target:
            cfg['num_epochs'] = compute_epochs(
                train_tf.num_triples, cfg['batch_size'], _target
            )
            logger.info(
                'target_steps=%d → num_epochs=%d (train triples=%d, batch=%d)',
                _target,
                cfg['num_epochs'],
                train_tf.num_triples,
                cfg['batch_size'],
            )
        else:
            raise ValueError("cfg must contain 'num_epochs' or 'target_steps'")

    # ── 2. Resolve device ─────────────────────────────────────────────────────
    if device == 'auto':
        if torch.cuda.is_available():
            device = 'cuda'
        elif torch.backends.mps.is_available():
            device = 'mps'
        else:
            device = 'cpu'

    if device == 'mps' and cfg.get('model') in _MPS_UNSUPPORTED_MODELS:
        logger.warning(
            '%s is not supported on MPS (Apple Silicon GPU). Falling back to CPU.',
            cfg['model'],
        )
        device = 'cpu'

    logger.info('Device: %s', device)

    # ── 3. Run PyKEEN pipeline ────────────────────────────────────────────────
    logger.info('=' * 70)
    logger.info(
        'Step 2/3 - Training %s  (epochs=%d, dim=%d, negs=%d)',
        cfg['model'],
        cfg['num_epochs'],
        cfg['embedding_dim'],
        cfg['num_negs_per_pos'],
    )
    logger.info('=' * 70)

    checkpoint_dir = out_dir / 'checkpoints'
    checkpoint_dir.mkdir(exist_ok=True)

    _es_patience = cfg.get('es_patience', 0)
    if _es_patience > 0:
        _es_frequency = cfg.get('es_frequency', 10)
        _es_delta = cfg.get('es_delta', 2e-3)
        stopper = 'early'
        stopper_kwargs = dict(
            metric='mrr',
            patience=_es_patience,
            frequency=_es_frequency,
            relative_delta=_es_delta,
            larger_is_better=True,
        )
        logger.info(
            'Early stopping - metric=mrr  patience=%d  frequency=%d  delta=%.4f',
            _es_patience,
            _es_frequency,
            _es_delta,
        )
    else:
        stopper = 'nop'
        stopper_kwargs = {}

    result = pipeline(
        training=train_tf,
        validation=val_tf,
        testing=test_tf,
        model=cfg['model'],
        model_kwargs=dict(embedding_dim=cfg['embedding_dim']),
        **({'loss': _MODEL_LOSS[cfg['model']]} if cfg['model'] in _MODEL_LOSS else {}),
        training_loop='sLCWA',
        negative_sampler='basic',
        negative_sampler_kwargs=dict(num_negs_per_pos=cfg['num_negs_per_pos']),
        optimizer='Adam',
        optimizer_kwargs=dict(lr=cfg['lr']),
        training_kwargs=dict(
            num_epochs=cfg['num_epochs'],
            batch_size=cfg['batch_size'],
            checkpoint_name=f'{cfg["model"].lower()}.pt',
            checkpoint_directory=str(checkpoint_dir),
            checkpoint_frequency=cfg['checkpoint_frequency'],
        ),
        # filter_known=True: when ranking (h,r,t), remove all other known true
        # tails for (h,r,?) from the candidate set. Standard KGE evaluation
        # protocol - intentional design choice, not a bug.
        evaluator='RankBasedEvaluator',
        evaluator_kwargs=dict(filtered=True),
        evaluation_kwargs=dict(batch_size=cfg['eval_batch_size']),
        stopper=stopper,
        stopper_kwargs=stopper_kwargs,
        device=device,
        random_seed=cfg['random_seed'],
        result_tracker='csv',
        result_tracker_kwargs=dict(path=str(out_dir / 'training_log.csv')),
    )

    # ── 4. Save artifacts ─────────────────────────────────────────────────────
    logger.info('=' * 60)
    logger.info('Step 3/3 - Saving model and results')
    result.save_to_directory(str(out_dir))
    test_tf.to_path_binary(out_dir / 'test_triples')
    logger.info('Saved PyKEEN pipeline result to: %s', out_dir)

    gp_metrics = _evaluate_gp_ranking(
        model=result.model,
        train_tf=train_tf,
        val_tf=val_tf,
        test_tf=test_tf,
        target_relation=target_relation,
        device=device,
        eval_batch_size=cfg['eval_batch_size'],
    )

    # ── 5. Build and save summary ─────────────────────────────────────────────
    summary = {
        'model': cfg['model'],
        'species': species,
        'target_relation': target_relation,
        'device': device,
        'include_relations': sorted(include_relations) if include_relations else None,
        'pre_split': pre_split,
        'graph': {
            'num_entities': _num_entities,
            'num_relations': _num_relations,
            'num_triples_total': _num_triples,
            'num_triples_train': train_tf.num_triples,
            'num_triples_val': val_tf.num_triples,
            'num_triples_test': test_tf.num_triples,
        },
        'config': cfg,
        'metrics': {
            'hits_at_1': _safe_metric(result, 'hits@1'),
            'hits_at_3': _safe_metric(result, 'hits@3'),
            'hits_at_10': _safe_metric(result, 'hits@10'),
            'mrr': _safe_metric(result, 'inverse_harmonic_mean_rank'),
            'mr': _safe_metric(result, 'mean_rank'),
        },
        'gp_metrics': gp_metrics,
    }

    summary_path = out_dir / 'summary.json'
    summary_path.write_text(json.dumps(summary, indent=2, default=str))

    logger.info('-' * 70)
    logger.info('Evaluation results (filtered, realistic):')
    for k, v in summary['metrics'].items():
        logger.info('  %s: %s', k, f'{v:.4f}' if isinstance(v, float) else v)
    if gp_metrics:
        logger.info('-' * 70)
        logger.info(
            'G→P 24-way tail ranking (%d candidates):',
            gp_metrics.get('gp_n_candidates', '?'),
        )
        for k in ('gp_mrr', 'gp_hits_at_1', 'gp_hits_at_3', 'gp_hits_at_10', 'gp_mr'):
            v = gp_metrics.get(k)
            logger.info('  %s: %s', k, f'{v:.4f}' if isinstance(v, float) else v)
    logger.info('Summary → %s', summary_path)

    return summary


# ── Helpers ───────────────────────────────────────────────────────────────────


def _evaluate_gp_ranking(
    model,
    train_tf: TriplesFactory,
    val_tf: TriplesFactory,
    test_tf: TriplesFactory,
    target_relation: str,
    device: str,
    eval_batch_size: int,
) -> dict:
    """
    24-way tail-restricted ranking evaluation on gene→phenotype triples only.

    Fixes the candidate set to the MP top-terms seen during training (~24 entities),
    making the task a constant-difficulty classification problem whose metrics are
    directly comparable across ablation conditions regardless of graph size.

    MRR/Hits@K with the unrestricted all-entity pool are NOT comparable across
    conditions because removing a layer (e.g. PPI) shrinks the entity pool and
    artificially inflates ranking metrics.
    """
    rel_id_test = test_tf.relation_to_id.get(target_relation)
    if rel_id_test is None:
        logger.warning(
            "G→P relation '%s' absent from test set - skipping 24-way eval",
            target_relation,
        )
        return {}

    gp_mask = test_tf.mapped_triples[:, 1] == rel_id_test
    gp_test = test_tf.mapped_triples[gp_mask]
    if len(gp_test) == 0:
        logger.warning('No G→P triples in test fold - skipping 24-way eval')
        return {}

    rel_id_train = train_tf.relation_to_id.get(target_relation)
    if rel_id_train is not None:
        train_gp_mask = train_tf.mapped_triples[:, 1] == rel_id_train
        mp_ids = torch.unique(train_tf.mapped_triples[train_gp_mask][:, 2])
    else:
        logger.warning(
            'G→P relation absent from training set - using test MP terms as candidates'
        )
        mp_ids = torch.unique(gp_test[:, 2])

    n_candidates = int(mp_ids.shape[0])
    n_test = int(gp_test.shape[0])
    logger.info(
        'G→P 24-way eval - %d test triples | %d MP candidates', n_test, n_candidates
    )

    # PyKEEN requires all entity IDs in mapped_triples (heads + tails) to be
    # present in restrict_entities_to, so include gene head IDs alongside the
    # 24 MP tail candidates.
    gene_ids_in_test = torch.unique(gp_test[:, 0])
    restrict_ids = torch.cat([mp_ids, gene_ids_in_test]).unique()

    evaluator = RankBasedEvaluator(filtered=True)
    try:
        eval_result = evaluator.evaluate(
            model=model,
            mapped_triples=gp_test,
            additional_filter_triples=[
                train_tf.mapped_triples,
                val_tf.mapped_triples,
            ],
            restrict_entities_to=restrict_ids,
            batch_size=eval_batch_size,
            device=device,
        )
    except Exception:
        logger.exception('24-way G→P evaluation failed')
        return {}

    def _get(primary: str, fallback: str | None = None) -> float | None:
        """
        Extract a metric by name, trying a fallback name on failure.
        """
        try:
            return float(eval_result.get_metric(primary))
        except Exception:
            pass
        if fallback:
            try:
                return float(eval_result.get_metric(fallback))
            except Exception:
                pass
        return None

    return {
        'gp_mrr': _get(
            'tail.realistic.inverse_harmonic_mean_rank', 'inverse_harmonic_mean_rank'
        ),
        'gp_hits_at_1': _get('tail.realistic.hits@1', 'hits@1'),
        'gp_hits_at_3': _get('tail.realistic.hits@3', 'hits@3'),
        'gp_hits_at_10': _get('tail.realistic.hits@10', 'hits@10'),
        'gp_mr': _get('tail.realistic.mean_rank', 'mean_rank'),
        'gp_n_candidates': n_candidates,
        'gp_n_test': n_test,
    }


def _safe_metric(result, key: str) -> float | None:
    """
    Extract a metric from a PyKEEN PipelineResult, returning None on failure.
    """
    try:
        return float(result.get_metric(key))
    except Exception:
        return None


def _report_relation_split(
    train_tf: TriplesFactory,
    val_tf: TriplesFactory,
    test_tf: TriplesFactory,
    target_relation: str,
) -> None:
    """
    Log the count of target-relation triples in each split.
    """
    rel_id = train_tf.relation_to_id.get(target_relation)
    if rel_id is None:
        logger.warning(
            "Target relation '%s' not found in training triples; split report skipped",
            target_relation,
        )
        return

    def _count(tf: TriplesFactory) -> int:
        """
        Count triples with the target relation in a TriplesFactory.
        """
        return int((tf.mapped_triples[:, 1] == rel_id).sum().item())

    logger.info(
        'Gene→Phenotype triples - train: %d  val: %d  test: %d',
        _count(train_tf),
        _count(val_tf),
        _count(test_tf),
    )


# ── CLI entry point ───────────────────────────────────────────────────────────

_FAST_PRESET = {
    'embedding_dim': 64,
    'num_epochs': 30,
    'num_negs_per_pos': 32,
    'batch_size': 1024,
}


def _parse_args() -> argparse.Namespace:
    """
    Parse and return CLI arguments for the KGE train entry point.
    """
    p = argparse.ArgumentParser(
        description='Train a KGE model on the knowledge graph (gene→phenotype focus).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        '--data-dir',
        required=True,
        help='BioCypher output dir with *-header.csv and *-part000.csv files.',
    )
    p.add_argument(
        '--out-dir',
        default=None,
        help='Directory for model artifacts. Defaults to pykeen_out/<model_lower>/.',
    )
    p.add_argument(
        '--species',
        default='human',
        choices=['human', 'mouse'],
        help='Species used to select the target gene→phenotype relation.',
    )
    p.add_argument(
        '--target-relation',
        default=None,
        help='Explicit target relation. Overrides --species.',
    )
    p.add_argument('--epochs', type=int, default=DEFAULT_CONFIG.get('num_epochs'))
    p.add_argument('--dim', type=int, default=DEFAULT_CONFIG['embedding_dim'])
    p.add_argument('--lr', type=float, default=DEFAULT_CONFIG['lr'])
    p.add_argument('--batch', type=int, default=DEFAULT_CONFIG['batch_size'])
    p.add_argument('--negs', type=int, default=DEFAULT_CONFIG['num_negs_per_pos'])
    p.add_argument(
        '--model', default='RotatE', choices=['RotatE', 'ComplEx', 'DistMult', 'TransE']
    )
    p.add_argument('--device', default='auto', choices=['auto', 'cpu', 'cuda', 'mps'])
    p.add_argument('--seed', type=int, default=DEFAULT_CONFIG['random_seed'])
    p.add_argument(
        '--fast',
        action='store_true',
        help='Smoke-test preset: dim=64, epochs=30, negs=32, batch=1024.',
    )
    p.add_argument('--cache-dir', default=None)
    return p.parse_args()


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s  %(levelname)-8s  %(name)s  %(message)s',
        datefmt='%H:%M:%S',
    )
    args = _parse_args()

    if args.out_dir is None:
        args.out_dir = f'pykeen_out/{args.model.lower()}'

    override = (
        _FAST_PRESET.copy()
        if args.fast
        else {
            'num_epochs': args.epochs,
            'embedding_dim': args.dim,
            'lr': args.lr,
            'batch_size': args.batch,
            'num_negs_per_pos': args.negs,
            'random_seed': args.seed,
        }
    )
    override['model'] = args.model

    if args.model in _COMPLEX_SPACE_MODELS:
        dim = override.get('embedding_dim', args.dim)
        if dim % 2 != 0:
            raise ValueError(
                f'embedding_dim must be even for {args.model} (complex space, got {dim})'
            )

    run_training(
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        config=override,
        species=args.species,
        target_relation=args.target_relation,
        device=args.device,
        cache_dir=args.cache_dir,
    )
