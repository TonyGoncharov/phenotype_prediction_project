"""
src/machine_learning/common/splits.py - reproducible data splits shared by both ML tracks.

Provides gene-level splits for the R-GCN classifier track. PyKEEN TriplesFactory
splitting lives in kge/train.py because it depends on PyKEEN internals.
"""

from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger(__name__)


def make_gene_split(
    all_genes: set[str],
    train_frac: float = 0.8,
    seed: int = 42,
) -> tuple[list[str], list[str]]:
    """
    Return a reproducible (train_genes, test_genes) gene-level split.

    Genes are sorted before shuffling so the split is deterministic across
    Python versions and platforms.

    Args:
        all_genes: Complete set of gene IDs to split.
        train_frac: Fraction of genes assigned to training.
        seed: NumPy RNG seed.

    Returns:
        (train_genes, test_genes) - disjoint lists covering all_genes.
    """
    genes = sorted(all_genes)
    rng = np.random.default_rng(seed)
    idx = np.arange(len(genes))
    rng.shuffle(idx)
    n_train = int(len(idx) * train_frac)

    train_genes = [genes[i] for i in idx[:n_train]]
    test_genes = [genes[i] for i in idx[n_train:]]

    logger.info(
        'Gene split (seed=%d): %d train / %d test (total %d)',
        seed,
        len(train_genes),
        len(test_genes),
        len(genes),
    )
    return train_genes, test_genes
