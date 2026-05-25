"""
src/machine_learning/rgcn/train.py - R-GCN encoder training via link reconstruction.

Trains an RGCNEncoder on context edges (GO, expression, PPI, gene-protein)
using a self-supervised link-reconstruction objective with binary cross-entropy
loss. Phenotype edges are never seen during training.
"""

from __future__ import annotations

import logging

import torch
import torch.nn.functional as F
from src.machine_learning.rgcn.config import DEFAULT_CONFIG
from src.machine_learning.rgcn.model import RGCNEncoder

logger = logging.getLogger(__name__)


def train_rgcn(
    edge_index: torch.Tensor,
    edge_type: torch.Tensor,
    num_entities: int,
    num_relations: int,
    dim: int,
    epochs: int,
    device: torch.device | str,
    lr: float = DEFAULT_CONFIG['lr'],
    batch_size: int = DEFAULT_CONFIG['batch_size'],
    num_negs_per_pos: int = DEFAULT_CONFIG['num_negs_per_pos'],
    log_every: int = 20,
) -> RGCNEncoder:
    """
    Train an RGCNEncoder via link reconstruction on context edges.

    Args:
        edge_index: Bidirectional edge tensor, shape (2, 2 * n_edges).
        edge_type: Relation type indices, shape (2 * n_edges,).
        num_entities: Total number of entities.
        num_relations: Total number of relation types (including inverses).
        dim: Embedding dimensionality.
        epochs: Number of full training epochs.
        device: PyTorch device.
        lr: Adam learning rate.
        batch_size: Number of positive edges per gradient step.
        num_negs_per_pos: Negative samples per positive edge.
        log_every: Log training loss every this many epochs.

    Returns:
        Trained RGCNEncoder in eval mode.
    """
    device = torch.device(device)
    model = RGCNEncoder(num_entities, num_relations, dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=0.5,
        patience=10,
        min_lr=1e-5,
    )

    ei = edge_index.to(device)
    et = edge_type.to(device)
    pos_h = ei[0]
    pos_t = ei[1]

    logger.info('Training R-GCN for %d epochs, dim=%d …', epochs, dim)
    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad()

        z = model(ei, et)

        idx = torch.randperm(len(pos_h), device=device)[:batch_size]
        ph, pt = pos_h[idx], pos_t[idx]

        pos_scores = model.decode(z, ph, pt)

        neg_t = torch.randint(
            0, num_entities, (len(ph) * num_negs_per_pos,), device=device
        )
        neg_h = ph.repeat_interleave(num_negs_per_pos)
        neg_scores = model.decode(z, neg_h, neg_t)

        loss = F.binary_cross_entropy_with_logits(
            pos_scores, torch.ones_like(pos_scores)
        ) + F.binary_cross_entropy_with_logits(neg_scores, torch.zeros_like(neg_scores))
        loss.backward()
        optimizer.step()
        scheduler.step(loss.item())

        if epoch % log_every == 0:
            logger.info('  Epoch %3d/%d  loss=%.4f', epoch, epochs, loss.item())

    model.eval()
    return model
