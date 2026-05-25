"""
src/machine_learning/rgcn/evaluate.py - evaluation utilities for the R-GCN classification track.

Provides delta-metric computation for ablation comparisons. The primary metric
for R-GCN is mean AUPRC across phenotype classes. AUROC and per-class metrics
are also available via common/metrics.py.

Metric applicability note
-------------------------
MRR and Hits@K are KGE ranking metrics and do NOT apply to the R-GCN track.
R-GCN produces gene embeddings fed to a classifier; the output is class
probabilities, not triple scores, so ranking-based metrics are not meaningful.
"""

from __future__ import annotations

import logging

import pandas as pd
from src.machine_learning.common.metrics import compute_delta_metrics

logger = logging.getLogger(__name__)


def compute_ablation_deltas(
    comparison_df: pd.DataFrame,
    baseline_condition: str = 'all',
) -> pd.DataFrame:
    """
    Augment an ablation comparison DataFrame with delta columns.

    Delta = metric(all) − metric(condition). A positive delta means the
    full-graph condition outperformed, i.e. removing that layer hurt performance.

    Args:
        comparison_df: DataFrame with one row per condition and metric columns.
        baseline_condition: The condition name to treat as baseline.

    Returns:
        Input DataFrame with additional delta_* columns appended.
    """
    if baseline_condition not in comparison_df['condition'].values:
        logger.warning(
            "Baseline condition '%s' not found in comparison DataFrame - "
            'delta metrics will not be computed.',
            baseline_condition,
        )
        return comparison_df

    baseline_row = (
        comparison_df[comparison_df['condition'] == baseline_condition]
        .iloc[0]
        .to_dict()
    )

    delta_rows = []
    for _, row in comparison_df.iterrows():
        deltas = compute_delta_metrics(baseline_row, row.to_dict())
        delta_rows.append(deltas)

    delta_df = pd.DataFrame(delta_rows)
    return pd.concat([comparison_df.reset_index(drop=True), delta_df], axis=1)


def leakage_check(
    context_df: pd.DataFrame,
    pheno_df: pd.DataFrame,
    pheno_rel: str,
) -> None:
    """
    Assert that no phenotype edges appear in the context graph.

    Raises a RuntimeError if phenotype relation strings are found in context_df.
    Logs a warning for each violation found.

    Args:
        context_df: Context graph DataFrame (head, relation, tail).
        pheno_df: Phenotype edges DataFrame (head, relation, tail).
        pheno_rel: Phenotype relation name to check for.
    """
    if pheno_rel in context_df['relation'].values:
        n = (context_df['relation'] == pheno_rel).sum()
        raise RuntimeError(
            f'LEAKAGE DETECTED: {n} {pheno_rel!r} edges found in the context graph. '
            'These edges must never be in the R-GCN message-passing graph.'
        )

    gene_pheno_pairs = set(zip(pheno_df['head'], pheno_df['tail'], strict=True))
    n_leaked = 0
    for _, row in context_df.iterrows():
        if (row['head'], row['tail']) in gene_pheno_pairs:
            n_leaked += 1
    if n_leaked:
        logger.warning(
            'LEAKAGE WARNING: %d context edges share (head, tail) pairs with '
            'phenotype annotations. Check that indirect leakage paths are absent.',
            n_leaked,
        )
