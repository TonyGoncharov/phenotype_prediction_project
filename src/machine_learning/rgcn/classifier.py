"""
src/machine_learning/rgcn/classifier.py - multi-label phenotype classifier on R-GCN embeddings.

Trains a OneVsRest LogisticRegression on gene embeddings produced by the R-GCN
encoder. The classifier is trained on a fixed 80/20 gene split that is shared
across all ablation conditions, making AUPRC values directly comparable.

MultiLabelBinarizer is always fitted on the complete gene set before the split
so the class list (MP top terms) is identical in every condition.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import MultiLabelBinarizer
from src.machine_learning.rgcn.config import DEFAULT_CONFIG
from src.machine_learning.rgcn.embeddings import embed_genes

logger = logging.getLogger(__name__)


def build_label_matrix(
    gene_list: list[str],
    gene_to_pheno: dict[str, list[str]],
    mlb: MultiLabelBinarizer,
) -> np.ndarray:
    """
    Build a binary label matrix for *gene_list* using a fitted MultiLabelBinarizer.

    Args:
        gene_list: Ordered list of gene IDs.
        gene_to_pheno: Mapping of gene ID → list of MP term IDs.
        mlb: Already-fitted MultiLabelBinarizer (determines class order).

    Returns:
        Binary array of shape (len(gene_list), n_classes).
    """
    return mlb.transform([gene_to_pheno.get(g, []) for g in gene_list])


def train_classifier(
    X_train: np.ndarray,
    Y_train: np.ndarray,
    seed: int = DEFAULT_CONFIG['seed'],
    n_jobs: int = -1,
    max_iter: int = DEFAULT_CONFIG['classifier_max_iter'],
    C: float = DEFAULT_CONFIG['classifier_C'],
) -> OneVsRestClassifier:
    """
    Train a OneVsRest LogisticRegression on gene embeddings.

    Args:
        X_train: Gene embedding matrix, shape (n_train, dim).
        Y_train: Binary multi-label matrix, shape (n_train, n_classes).
        seed: Random seed for LogisticRegression.
        n_jobs: Parallelism. Use 1 when running inside a worker process to
            avoid double-forking; use -1 in sequential mode.
        max_iter: Maximum LogisticRegression iterations.
        C: Inverse regularisation strength for LogisticRegression.

    Returns:
        Fitted OneVsRestClassifier.
    """
    clf = OneVsRestClassifier(
        LogisticRegression(max_iter=max_iter, random_state=seed, C=C),
        n_jobs=n_jobs,
    )
    clf.fit(X_train, Y_train)
    return clf


def evaluate_classifier(
    clf: OneVsRestClassifier,
    X_test: np.ndarray,
    Y_test: np.ndarray,
    class_labels: list[str],
) -> tuple[float, pd.DataFrame]:
    """
    Evaluate a fitted classifier and return mean AUPRC + per-class metrics.

    Args:
        clf: Fitted OneVsRestClassifier.
        X_test: Gene embedding matrix, shape (n_test, dim).
        Y_test: Binary ground-truth matrix, shape (n_test, n_classes).
        class_labels: MP term IDs, length n_classes.

    Returns:
        (mean_auprc, per_class_df) - per_class_df has columns
        [mp_term_id, ap, n_pos_test, n_test, baseline]. Sorted by ap descending.
        Classes with zero test positives are omitted.
    """
    Y_prob = clf.predict_proba(X_test)

    rows = []
    for i, mp_term in enumerate(class_labels):
        y_true = Y_test[:, i]
        y_score = Y_prob[:, i]
        if y_true.sum() == 0:
            continue
        rows.append(
            {
                'mp_term_id': mp_term,
                'ap': round(float(average_precision_score(y_true, y_score)), 4),
                'n_pos_test': int(y_true.sum()),
                'n_test': int(len(y_true)),
                'baseline': round(float(y_true.mean()), 4),
            }
        )

    per_class = (
        pd.DataFrame(rows).sort_values('ap', ascending=False).reset_index(drop=True)
    )
    mean_auprc = float(per_class['ap'].mean()) if len(per_class) > 0 else float('nan')
    return mean_auprc, per_class


def run_classification(
    embeddings: np.ndarray,
    ent2id: dict[str, int],
    pheno_df: pd.DataFrame,
    all_genes: set[str],
    train_genes: list[str],
    test_genes: list[str],
    mlb: MultiLabelBinarizer,
    dim: int,
    seed: int = DEFAULT_CONFIG['seed'],
    n_jobs: int = -1,
    classifier_max_iter: int = DEFAULT_CONFIG['classifier_max_iter'],
    classifier_C: float = DEFAULT_CONFIG['classifier_C'],
) -> tuple[float, pd.DataFrame]:
    """
    Full classification pipeline: embed genes, train, evaluate.

    Genes absent from ent2id receive zero embeddings (they are not dropped).
    This is the canonical entry point used by the ablation worker.

    Args:
        embeddings: Entity embedding array from get_embeddings(), shape (num_entities, dim).
        ent2id: Entity string → integer index.
        pheno_df: DataFrame with columns [head, relation, tail] for phenotype edges.
        all_genes: All phenotype-annotated genes (determines the universe).
        train_genes: Ordered list of training gene IDs.
        test_genes: Ordered list of test gene IDs.
        mlb: Fitted MultiLabelBinarizer (class list fixed across conditions).
        dim: Embedding dimensionality (for zero-vector construction).
        seed: Random seed for classifier.
        n_jobs: Parallelism for classifier.
        classifier_max_iter: Maximum LogisticRegression iterations.
        classifier_C: Inverse regularisation strength for LogisticRegression.

    Returns:
        (mean_auprc, per_class_df)
    """
    gene_to_pheno: dict[str, list[str]] = {}
    for _, row in pheno_df.iterrows():
        gene_to_pheno.setdefault(row['head'], []).append(row['tail'])

    n_missing_train = sum(1 for g in train_genes if g not in ent2id)
    n_missing_test = sum(1 for g in test_genes if g not in ent2id)
    if n_missing_train or n_missing_test:
        logger.info(
            'Zero embeddings assigned: %d train genes, %d test genes not in context graph',
            n_missing_train,
            n_missing_test,
        )

    X_train = embed_genes(embeddings, ent2id, train_genes, dim)
    X_test = embed_genes(embeddings, ent2id, test_genes, dim)
    Y_train = build_label_matrix(train_genes, gene_to_pheno, mlb)
    Y_test = build_label_matrix(test_genes, gene_to_pheno, mlb)

    logger.info(
        'Classifier - train: %d genes | test: %d genes | classes: %d',
        len(train_genes),
        len(test_genes),
        len(mlb.classes_),
    )

    clf = train_classifier(
        X_train,
        Y_train,
        seed=seed,
        n_jobs=n_jobs,
        max_iter=classifier_max_iter,
        C=classifier_C,
    )
    return evaluate_classifier(clf, X_test, Y_test, list(mlb.classes_))
