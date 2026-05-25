"""
src/machine_learning/rgcn/ablation.py - layer ablation study for the R-GCN track.

Trains one R-GCN encoder per condition (removing one context layer at a time),
then evaluates AUPRC with a fixed gene split and fixed class list.
"""

from __future__ import annotations

import logging
import multiprocessing
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from src.machine_learning.common.data_loading import load_triples, resolve_biocypher_dir
from src.machine_learning.common.results import (
    comparison_dataframe,
    save_comparison_csv,
)
from src.machine_learning.common.splits import make_gene_split
from src.machine_learning.rgcn.config import DEFAULT_CONFIG, target_relation_for_species
from src.machine_learning.rgcn.graph_builder import (
    CONTEXT_RELS_HUMAN,
    CONTEXT_RELS_MOUSE,
)

logger = logging.getLogger(__name__)

# ── Thread environment variables - must be set before torch import ─────────────

_THREAD_ENV_VARS = [
    'OMP_NUM_THREADS',
    'MKL_NUM_THREADS',
    'OPENBLAS_NUM_THREADS',
    'NUMEXPR_NUM_THREADS',
    'VECLIB_MAXIMUM_THREADS',
    'PYTORCH_NUM_THREADS',
]


def _set_thread_env(n_threads: int) -> None:
    """
    ProcessPoolExecutor initializer: set all thread caps before torch import.
    """
    val = str(n_threads)
    for var in _THREAD_ENV_VARS:
        os.environ[var] = val


# ── Ablation conditions ───────────────────────────────────────────────────────

_H_ENCODES = {'HumanGeneEncodesProtein'}
_H_GO = {'HumanGeneHasGoTerm'}
_H_EXPR = {'HumanGeneExpressedInTissue', 'TissueMappedToUberon'}
_H_PPI = {'HumanProteinInteractsWith'}

_M_ENCODES = {'MouseGeneEncodesProtein'}
_M_GO = {'MouseGeneHasGoTerm'}
_M_EXPR = {'MouseGeneExpressedInTissue'}
_M_PPI = {'MouseProteinInteractsWith'}

HUMAN_CONDITIONS: dict[str, set[str] | None] = {
    'all': None,
    'no_go': _H_ENCODES | _H_EXPR | _H_PPI,
    'no_ppi': _H_ENCODES | _H_GO | _H_EXPR,
    'no_expr': _H_ENCODES | _H_GO | _H_PPI,
    'encodes_only': _H_ENCODES,
}

MOUSE_CONDITIONS: dict[str, set[str] | None] = {
    'all': None,
    'no_go': _M_ENCODES | _M_EXPR | _M_PPI,
    'no_ppi': _M_ENCODES | _M_GO | _M_EXPR,
    'no_expr': _M_ENCODES | _M_GO | _M_PPI,
    'encodes_only': _M_ENCODES,
}

_SPECIES_CONDITIONS = {'human': HUMAN_CONDITIONS, 'mouse': MOUSE_CONDITIONS}
_SPECIES_ALL_CONTEXT = {'human': CONTEXT_RELS_HUMAN, 'mouse': CONTEXT_RELS_MOUSE}

_FAST_PRESET = {'rgcn_epochs': 20, 'dim': 64}


# ── Single-condition worker (runs in subprocess) ──────────────────────────────


def run_condition_worker(
    cond_name: str,
    include_rels: set[str] | None,
    all_df_rows: list,
    gene_to_pheno: dict[str, list[str]],
    train_genes: list[str],
    test_genes: list[str],
    mlb_classes: list[str],
    all_genes: set[str],
    out_dir: str,
    species: str,
    target_relation: str,
    rgcn_epochs: int,
    dim: int,
    seed: int,
    log_file: str,
    n_jobs: int = -1,
    lr: float = DEFAULT_CONFIG['lr'],
    batch_size: int = DEFAULT_CONFIG['batch_size'],
    num_negs_per_pos: int = DEFAULT_CONFIG['num_negs_per_pos'],
    classifier_max_iter: int = DEFAULT_CONFIG['classifier_max_iter'],
    classifier_C: float = DEFAULT_CONFIG['classifier_C'],
) -> dict:
    """
    Train R-GCN for one ablation condition and return a result dict.

    This function is designed to run inside a spawned worker process. Heavy
    imports (torch, sklearn) are deferred to after this function is called,
    because the spawn initializer (_set_thread_env) must run first to configure
    OpenMP thread counts before PyTorch initialises its thread pool.

    Args:
        cond_name: Condition identifier string.
        include_rels: Active context relation set (None → species-default full set).
        all_df_rows: Serialised triples as list of [head, relation, tail] lists.
        gene_to_pheno: Gene ID → list of MP term IDs.
        train_genes: Fixed training gene IDs.
        test_genes: Fixed test gene IDs.
        mlb_classes: Sorted MP class list (same across all conditions).
        all_genes: All phenotype-annotated genes.
        out_dir: Root output directory path string.
        species: 'human' or 'mouse'.
        target_relation: Gene→phenotype relation used as labels.
        rgcn_epochs: Number of R-GCN training epochs.
        dim: Embedding dimensionality.
        seed: Random seed.
        log_file: Per-condition log file path.
        n_jobs: Parallelism for classifier (-1 in sequential mode, 1 in parallel).
        lr: Adam learning rate for R-GCN training.
        batch_size: Positive edges sampled per R-GCN gradient step.
        num_negs_per_pos: Negative tails per positive edge.
        classifier_max_iter: Maximum LogisticRegression iterations.
        classifier_C: Inverse regularisation strength for LogisticRegression.

    Returns:
        Dict with condition metrics or {'condition': cond_name, 'error': str}.
    """
    import logging as _logging

    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    fh = _logging.FileHandler(log_file)
    fh.setLevel(_logging.INFO)
    fh.setFormatter(_logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    cond_log = _logging.getLogger(f'rgcn.{cond_name}')
    cond_log.setLevel(_logging.INFO)
    cond_log.addHandler(fh)
    cond_log.propagate = False

    t0 = time.time()
    cond_log.info("Condition '%s' started (pid=%d)", cond_name, os.getpid())
    cond_log.info('OMP_NUM_THREADS=%s', os.environ.get('OMP_NUM_THREADS', 'unset'))

    try:
        # Deferred imports - must happen AFTER thread env vars are set.
        import torch
        from sklearn.preprocessing import MultiLabelBinarizer
        from src.machine_learning.common.results import save_json
        from src.machine_learning.rgcn.classifier import run_classification
        from src.machine_learning.rgcn.embeddings import get_embeddings
        from src.machine_learning.rgcn.graph_builder import build_context_graph
        from src.machine_learning.rgcn.train import train_rgcn

        all_df = pd.DataFrame(all_df_rows, columns=['head', 'relation', 'tail'])
        pheno_df = all_df[all_df['relation'] == target_relation].copy()

        graph = build_context_graph(
            all_df,
            species=species,
            active_rels=include_rels,
            target_relation=target_relation,
        )

        n_missing = sum(1 for g in all_genes if g not in graph.ent2id)
        cond_log.info(
            'Graph: %d entities | %d rel types | %d edges | %d genes missing',
            graph.num_entities,
            graph.num_relations // 2,
            len(graph.context_df),
            n_missing,
        )

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model = train_rgcn(
            graph.edge_index,
            graph.edge_type,
            graph.num_entities,
            graph.num_relations,
            dim=dim,
            epochs=rgcn_epochs,
            device=device,
            lr=lr,
            batch_size=batch_size,
            num_negs_per_pos=num_negs_per_pos,
        )

        embeddings = get_embeddings(model, graph.edge_index, graph.edge_type, device)

        cond_path = Path(out_dir) / cond_name
        cond_path.mkdir(parents=True, exist_ok=True)
        torch.save(model.state_dict(), cond_path / 'rgcn_encoder.pt')

        mlb = MultiLabelBinarizer(classes=mlb_classes)
        mlb.fit([])

        mean_auprc, per_class = run_classification(
            embeddings=embeddings,
            ent2id=graph.ent2id,
            pheno_df=pheno_df,
            all_genes=all_genes,
            train_genes=train_genes,
            test_genes=test_genes,
            mlb=mlb,
            dim=dim,
            seed=seed,
            n_jobs=n_jobs,
            classifier_max_iter=classifier_max_iter,
            classifier_C=classifier_C,
        )

        per_class.to_csv(cond_path / 'auprc_per_phenotype.csv', index=False)

        mean_baseline = (
            float(per_class['baseline'].mean()) if len(per_class) > 0 else 0.0
        )
        improvement = mean_auprc / mean_baseline if mean_baseline > 0 else 0.0
        elapsed = time.time() - t0

        result = {
            'condition': cond_name,
            'species': species,
            'target_relation': target_relation,
            'relations': 'all'
            if include_rels is None
            else '|'.join(sorted(include_rels)),
            'mean_auprc': round(mean_auprc, 4),
            'mean_baseline': round(mean_baseline, 4),
            'improvement_x': round(improvement, 3),
            'n_entities': graph.num_entities,
            'n_context_edges': len(graph.context_df),
            'n_genes_in_graph': sum(1 for g in all_genes if g in graph.ent2id),
            'n_genes_missing': n_missing,
            'rgcn_epochs': rgcn_epochs,
            'rgcn_dim': dim,
            'lr': lr,
            'batch_size': batch_size,
            'num_negs_per_pos': num_negs_per_pos,
            'elapsed_s': round(elapsed, 1),
        }
        save_json(result, cond_path / 'summary.json')

        cond_log.info(
            "Condition '%s' done - AUPRC=%.4f  improvement=%.2fx  (%.0fs)",
            cond_name,
            mean_auprc,
            improvement,
            elapsed,
        )

    except Exception as exc:
        elapsed = time.time() - t0
        cond_log.exception("Condition '%s' FAILED", cond_name)
        result = {
            'condition': cond_name,
            'error': str(exc),
            'elapsed_s': round(elapsed, 1),
        }

    finally:
        cond_log.removeHandler(fh)
        fh.close()

    return result


# ── Shared data preparation ────────────────────────────────────────────────────


def _prepare_shared_data(
    data_dir: str,
    species: str,
    seed: int,
    train_frac: float,
    target_relation: str,
) -> tuple[list, dict, list[str], list[str], set[str], list[str]]:
    """
    Load triples and derive all shared inputs for the ablation loop.

    Returns data in pickle-safe forms (lists/dicts, not DataFrames or sklearn
    objects) so they transfer correctly through ProcessPoolExecutor.

    Returns:
        (all_df_rows, gene_to_pheno, train_genes, test_genes, all_genes, mlb_classes)
    """
    raw = load_triples(data_dir)
    all_df = pd.DataFrame(raw, columns=['head', 'relation', 'tail'])

    pheno_df = all_df[all_df['relation'] == target_relation].copy()
    all_genes: set[str] = set(pheno_df['head'].unique())

    logger.info('Phenotype-annotated genes: %d', len(all_genes))
    logger.info('Gene→Phenotype triples: %d', len(pheno_df))
    logger.info('Unique MP top-terms: %d', pheno_df['tail'].nunique())

    train_genes, test_genes = make_gene_split(
        all_genes,
        train_frac=train_frac,
        seed=seed,
    )
    logger.info(
        'Fixed gene split (seed=%d): %d train / %d test',
        seed,
        len(train_genes),
        len(test_genes),
    )

    gene_to_pheno: dict[str, list[str]] = {}
    for _, row in pheno_df.iterrows():
        gene_to_pheno.setdefault(row['head'], []).append(row['tail'])

    all_mp = sorted({mp for mps in gene_to_pheno.values() for mp in mps})
    logger.info('Phenotype classes (shared): %d', len(all_mp))

    return (
        all_df.values.tolist(),
        gene_to_pheno,
        train_genes,
        test_genes,
        all_genes,
        all_mp,
    )


# ── Parallel runner ────────────────────────────────────────────────────────────


def _run_parallel(
    conditions: dict[str, set[str] | None],
    shared: tuple,
    out_dir: str,
    species: str,
    target_relation: str,
    rgcn_epochs: int,
    dim: int,
    seed: int,
    max_workers: int,
    threads_per_worker: int,
    lr: float,
    batch_size: int,
    num_negs_per_pos: int,
    classifier_max_iter: int,
    classifier_C: float,
) -> list[dict]:
    """
    Run all conditions in parallel using ProcessPoolExecutor with spawn context.

    spawn is mandatory: it starts a fresh interpreter per worker so that
    _set_thread_env runs before any torch import. With fork, torch would already
    be imported and the OpenMP pool already initialised - setting env vars
    afterwards would have no effect.
    """
    all_df_rows, gene_to_pheno, train_genes, test_genes, all_genes, mlb_classes = shared

    logger.info(
        'Parallel mode: %d workers × %d threads each',
        max_workers,
        threads_per_worker,
    )

    futures = {}
    results = []

    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=multiprocessing.get_context('spawn'),
        initializer=_set_thread_env,
        initargs=(threads_per_worker,),
    ) as pool:
        for cond_name, include_rels in conditions.items():
            log_file = str(Path(out_dir) / cond_name / 'worker.log')
            future = pool.submit(
                run_condition_worker,
                cond_name=cond_name,
                include_rels=include_rels,
                all_df_rows=all_df_rows,
                gene_to_pheno=gene_to_pheno,
                train_genes=train_genes,
                test_genes=test_genes,
                mlb_classes=mlb_classes,
                all_genes=all_genes,
                out_dir=out_dir,
                species=species,
                target_relation=target_relation,
                rgcn_epochs=rgcn_epochs,
                dim=dim,
                seed=seed,
                log_file=log_file,
                n_jobs=1,
                lr=lr,
                batch_size=batch_size,
                num_negs_per_pos=num_negs_per_pos,
                classifier_max_iter=classifier_max_iter,
                classifier_C=classifier_C,
            )
            futures[future] = cond_name
            logger.info('Submitted: %s  (log → %s)', cond_name, log_file)

        for future in as_completed(futures):
            cond_name = futures[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as exc:
                logger.exception("Condition '%s' raised an exception", cond_name)
                results.append({'condition': cond_name, 'error': str(exc)})

    return results


# ── Sequential runner ─────────────────────────────────────────────────────────


def _run_sequential(
    conditions: dict[str, set[str] | None],
    shared: tuple,
    out_dir: str,
    species: str,
    target_relation: str,
    rgcn_epochs: int,
    dim: int,
    seed: int,
    lr: float,
    batch_size: int,
    num_negs_per_pos: int,
    classifier_max_iter: int,
    classifier_C: float,
) -> list[dict]:
    """
    Run each ablation condition sequentially in the current process.

    Args:
        conditions: Mapping of condition name → active relation set (None = all).
        shared: Tuple of pre-loaded shared data from _prepare_shared_data.
        out_dir: Root output directory; per-condition subdirs are created within.
        species: 'human' or 'mouse'.
        target_relation: Gene→phenotype relation name used as labels.
        rgcn_epochs: Number of R-GCN training epochs per condition.
        dim: R-GCN embedding dimension.
        seed: Random seed for training and classifier.
        lr: Learning rate for the R-GCN encoder.
        batch_size: Mini-batch size for link-reconstruction training.
        num_negs_per_pos: Negative samples per positive triple.
        classifier_max_iter: Maximum iterations for logistic regression.
        classifier_C: Regularisation strength for logistic regression.

    Returns:
        List of per-condition result dicts (one per condition).
    """
    all_df_rows, gene_to_pheno, train_genes, test_genes, all_genes, mlb_classes = shared
    results = []

    for cond_name, include_rels in conditions.items():
        rels_desc = 'all' if include_rels is None else '|'.join(sorted(include_rels))
        logger.info('=' * 70)
        logger.info('Condition: %s', cond_name)
        logger.info('Relations: %s', rels_desc)
        logger.info('=' * 70)

        log_file = str(Path(out_dir) / cond_name / 'worker.log')
        result = run_condition_worker(
            cond_name=cond_name,
            include_rels=include_rels,
            all_df_rows=all_df_rows,
            gene_to_pheno=gene_to_pheno,
            train_genes=train_genes,
            test_genes=test_genes,
            mlb_classes=mlb_classes,
            all_genes=all_genes,
            out_dir=out_dir,
            species=species,
            target_relation=target_relation,
            rgcn_epochs=rgcn_epochs,
            dim=dim,
            seed=seed,
            log_file=log_file,
            n_jobs=-1,
            lr=lr,
            batch_size=batch_size,
            num_negs_per_pos=num_negs_per_pos,
            classifier_max_iter=classifier_max_iter,
            classifier_C=classifier_C,
        )
        results.append(result)

    return results


# ── Main orchestrator ─────────────────────────────────────────────────────────


def run_ablation(
    data_dir: str | None = None,
    out_dir: str | None = None,
    conditions: dict[str, set[str] | None] | None = None,
    species: str = 'human',
    rgcn_epochs: int = DEFAULT_CONFIG['rgcn_epochs'],
    dim: int = DEFAULT_CONFIG['dim'],
    seed: int = DEFAULT_CONFIG['seed'],
    train_frac: float = DEFAULT_CONFIG['train_frac'],
    target_relation: str | None = None,
    lr: float = DEFAULT_CONFIG['lr'],
    batch_size: int = DEFAULT_CONFIG['batch_size'],
    num_negs_per_pos: int = DEFAULT_CONFIG['num_negs_per_pos'],
    classifier_max_iter: int = DEFAULT_CONFIG['classifier_max_iter'],
    classifier_C: float = DEFAULT_CONFIG['classifier_C'],
    parallel: bool = False,
    max_workers: int | None = None,
    threads_per_worker: int | None = None,
) -> pd.DataFrame:
    """
    Train and evaluate one R-GCN model per ablation condition.

    All conditions use the same fixed gene split, same fixed class list, and
    same negative sampling strategy (random tail corruption during encoding).
    AUPRC is the primary comparison metric.

    Args:
        data_dir: BioCypher output directory (None → auto-detect).
        out_dir: Root output directory. Defaults to results/ml/{species}/rgcn/.
        conditions: Condition name → relation set. None → species default.
        species: 'human' or 'mouse'.
        rgcn_epochs: Number of R-GCN training epochs per condition.
        dim: R-GCN embedding dimensionality.
        seed: Random seed for splits, sampling, and classifier.
        train_frac: Fraction of phenotype-annotated genes used for classifier training.
        target_relation: Explicit gene→phenotype relation. Defaults by species.
        lr: Adam learning rate for R-GCN training.
        batch_size: Positive edges sampled per R-GCN gradient step.
        num_negs_per_pos: Negative tails per positive edge.
        classifier_max_iter: Maximum LogisticRegression iterations.
        classifier_C: Inverse regularisation strength for LogisticRegression.
        parallel: Run conditions in parallel (spawn-based subprocesses).
        max_workers: Max parallel workers. Defaults to len(conditions).
        threads_per_worker: OMP_NUM_THREADS per worker. Defaults to cpu_count // n_workers.

    Returns:
        DataFrame with one row per condition, written to <out_dir>/comparison_summary.csv.
    """
    target_relation = target_relation_for_species(species, target_relation)
    resolved_dir = resolve_biocypher_dir(data_dir, species=species)
    conditions = conditions or _SPECIES_CONDITIONS.get(species, HUMAN_CONDITIONS)
    if out_dir is None:
        out_dir = f'results/ml/{species}/rgcn'
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    shared = _prepare_shared_data(
        str(resolved_dir),
        species=species,
        seed=seed,
        train_frac=train_frac,
        target_relation=target_relation,
    )

    if parallel:
        n_workers = max_workers or len(conditions)
        n_threads = threads_per_worker or max(1, (os.cpu_count() or 1) // n_workers)
        results = _run_parallel(
            conditions=conditions,
            shared=shared,
            out_dir=out_dir,
            species=species,
            target_relation=target_relation,
            rgcn_epochs=rgcn_epochs,
            dim=dim,
            seed=seed,
            max_workers=n_workers,
            threads_per_worker=n_threads,
            lr=lr,
            batch_size=batch_size,
            num_negs_per_pos=num_negs_per_pos,
            classifier_max_iter=classifier_max_iter,
            classifier_C=classifier_C,
        )
    else:
        results = _run_sequential(
            conditions=conditions,
            shared=shared,
            out_dir=out_dir,
            species=species,
            target_relation=target_relation,
            rgcn_epochs=rgcn_epochs,
            dim=dim,
            seed=seed,
            lr=lr,
            batch_size=batch_size,
            num_negs_per_pos=num_negs_per_pos,
            classifier_max_iter=classifier_max_iter,
            classifier_C=classifier_C,
        )

    df = comparison_dataframe(results)
    save_comparison_csv(df, out_dir)
    return df
