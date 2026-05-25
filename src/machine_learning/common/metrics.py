"""
src/machine_learning/common/metrics.py - metric computation shared across both ML tracks.

Both KGE (scored via predictor) and R-GCN (scored via classifier) ultimately
produce per-class AUPRC values that need to be compared across ablation conditions.
This module provides the shared aggregation and delta-metric logic.

Metric applicability:
  KGE track  : AUPRC (primary), AUROC, Precision@K, Recall@K, MRR, Hits@1/3/10
  R-GCN track : AUPRC (primary), AUROC, Precision@K, Recall@K
  MRR and Hits@K are ranking metrics specific to the KGE triple-scoring protocol.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

logger = logging.getLogger(__name__)

STANDARD_DELTA_METRICS = [
    'mean_auprc',
    'auroc',
    'mrr',
    'hits_at_1',
    'hits_at_3',
    'hits_at_10',
    'precision_at_k',
    'recall_at_k',
]


def per_class_auprc(
    y_true_matrix: np.ndarray,
    y_score_matrix: np.ndarray,
    class_labels: list[str],
    include_auroc: bool = True,
) -> pd.DataFrame:
    """
    Compute per-class average precision (and optionally AUROC) for multi-label output.

    Args:
        y_true_matrix: Binary ground-truth labels, shape (n_samples, n_classes).
        y_score_matrix: Predicted scores, shape (n_samples, n_classes).
        class_labels: Class name strings, length n_classes.
        include_auroc: Whether to compute AUROC alongside AUPRC.

    Returns:
        DataFrame with columns [mp_term_id, ap, auroc, n_pos_test, n_test, baseline].
        Sorted by ap descending. Classes with zero positives in y_true are omitted.
    """
    rows = []
    for i, label in enumerate(class_labels):
        y_true = y_true_matrix[:, i]
        y_score = y_score_matrix[:, i]
        n_pos = int(y_true.sum())
        if n_pos == 0:
            continue

        ap = float(average_precision_score(y_true, y_score))

        auroc = float('nan')
        if include_auroc:
            try:
                auroc = float(roc_auc_score(y_true, y_score))
            except ValueError as exc:
                logger.warning('AUROC undefined for class %s: %s', label, exc)

        rows.append(
            {
                'mp_term_id': label,
                'ap': round(ap, 4),
                'auroc': round(auroc, 4),
                'n_pos_test': n_pos,
                'n_test': int(len(y_true)),
                'baseline': round(float(y_true.mean()), 4),
            }
        )

    return pd.DataFrame(rows).sort_values('ap', ascending=False).reset_index(drop=True)


def compute_delta_metrics(
    baseline: dict,
    condition: dict,
    metric_keys: list[str] | None = None,
    suffix: str = '',
) -> dict:
    """
    Return Δmetric = metric(baseline) − metric(condition) for each key.

    A positive delta means the baseline ('all' condition) outperformed the
    ablated condition, so that layer contributed positively.

    Args:
        baseline: Metric dict for the full-graph ('all') condition.
        condition: Metric dict for the ablated condition.
        metric_keys: Keys to compute deltas for. Defaults to the standard set.
        suffix: Optional suffix for delta column names, e.g. '_vs_all'.

    Returns:
        Dict of {f'delta_{key}{suffix}': float | None}.
    """
    if metric_keys is None:
        metric_keys = STANDARD_DELTA_METRICS

    deltas: dict = {}
    for key in metric_keys:
        b = baseline.get(key)
        c = condition.get(key)
        delta_key = f'delta_{key}{suffix}'
        if b is not None and c is not None and not pd.isna(b) and not pd.isna(c):
            try:
                deltas[delta_key] = round(float(b) - float(c), 4)
            except (TypeError, ValueError):
                deltas[delta_key] = None
        else:
            deltas[delta_key] = None
    return deltas


def precision_at_k(y_true: np.ndarray, y_score: np.ndarray, k: int) -> float:
    """
    Precision at rank k for a single query.

    Args:
        y_true: Binary relevance labels.
        y_score: Predicted scores (higher = more relevant).
        k: Cutoff rank.

    Returns:
        Fraction of top-k retrieved items that are relevant.
    """
    top_k = np.argsort(y_score)[::-1][:k]
    return float(y_true[top_k].mean())


def recall_at_k(y_true: np.ndarray, y_score: np.ndarray, k: int) -> float:
    """
    Recall at rank k for a single query.

    Args:
        y_true: Binary relevance labels.
        y_score: Predicted scores (higher = more relevant).
        k: Cutoff rank.

    Returns:
        Fraction of all relevant items that appear in the top-k.
    """
    n_pos = int(y_true.sum())
    if n_pos == 0:
        return 0.0
    top_k = np.argsort(y_score)[::-1][:k]
    return float(y_true[top_k].sum()) / n_pos
