"""
src/machine_learning/rgcn/embeddings.py - gene embedding extraction, saving, and loading.

After R-GCN training, gene embeddings are extracted from the encoder and used
as input features for the downstream multi-label classifier.

Genes absent from the context graph (no context edges after ablation filtering)
receive zero embeddings - they are NOT dropped. This is an honest penalty for
missing information: a gene that the graph knows nothing about should not be
silently excluded from evaluation.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import torch
from src.machine_learning.rgcn.model import RGCNEncoder

logger = logging.getLogger(__name__)


def get_embeddings(
    model: RGCNEncoder,
    edge_index: torch.Tensor,
    edge_type: torch.Tensor,
    device: torch.device | str,
) -> np.ndarray:
    """
    Extract entity embedding matrix from a trained R-GCN encoder.

    Args:
        model: Trained RGCNEncoder (will be set to eval mode).
        edge_index: Bidirectional edge tensor, shape (2, 2 * n_edges).
        edge_type: Relation type indices, shape (2 * n_edges,).
        device: Device for inference.

    Returns:
        Float32 numpy array of shape (num_entities, hidden_dim).
    """
    device = torch.device(device)
    model.eval()
    with torch.no_grad():
        z = model(edge_index.to(device), edge_type.to(device))
    return z.cpu().numpy().astype(np.float32)


def embed_genes(
    embeddings: np.ndarray,
    ent2id: dict[str, int],
    gene_list: list[str],
    dim: int,
) -> np.ndarray:
    """
    Build a gene embedding matrix respecting the zero-embedding contract.

    Genes present in ent2id get their learned embedding. Genes absent from the
    context graph (e.g. after aggressive ablation filtering) receive a zero vector.

    Args:
        embeddings: Full entity embedding array, shape (num_entities, dim).
        ent2id: Entity string → integer index mapping.
        gene_list: Ordered list of gene IDs to embed.
        dim: Embedding dimensionality (used to construct zero vectors).

    Returns:
        Float32 array of shape (len(gene_list), dim).
    """
    rows = []
    for g in gene_list:
        rows.append(
            embeddings[ent2id[g]] if g in ent2id else np.zeros(dim, dtype=np.float32)
        )
    return np.array(rows, dtype=np.float32)


def save_encoder(
    model: RGCNEncoder,
    path: Path | str,
    ent2id: dict[str, int] | None = None,
    rel2id: dict[str, int] | None = None,
) -> None:
    """
    Save encoder state dict and optional entity/relation mappings.

    Args:
        model: Trained RGCNEncoder.
        path: Output file path (e.g. results/ml/human/rgcn/all/rgcn_encoder.pt).
        ent2id: Entity-to-index mapping to save alongside the weights.
        rel2id: Relation-to-index mapping to save alongside the weights.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict = {'state_dict': model.state_dict()}
    if ent2id is not None:
        payload['ent2id'] = ent2id
    if rel2id is not None:
        payload['rel2id'] = rel2id
    torch.save(payload, path)
    logger.info('Encoder saved → %s', path)


def load_encoder(
    path: Path | str,
    num_entities: int,
    num_relations: int,
    dim: int,
    device: torch.device | str = 'cpu',
) -> tuple[RGCNEncoder, dict[str, int] | None, dict[str, int] | None]:
    """
    Load an encoder from a checkpoint written by save_encoder().

    Args:
        path: Checkpoint file path.
        num_entities: Must match the checkpoint.
        num_relations: Must match the checkpoint.
        dim: Embedding dimension (must match the checkpoint).
        device: Device to load tensors onto.

    Returns:
        (model, ent2id, rel2id) - mappings are None if not stored in the checkpoint.
    """
    device = torch.device(device)
    payload = torch.load(path, map_location=device, weights_only=False)

    if isinstance(payload, dict) and 'state_dict' in payload:
        state_dict = payload['state_dict']
        ent2id = payload.get('ent2id')
        rel2id = payload.get('rel2id')
    else:
        # Legacy format: payload is the state_dict directly.
        state_dict = payload
        ent2id = None
        rel2id = None

    model = RGCNEncoder(num_entities, num_relations, dim).to(device)
    model.load_state_dict(state_dict)
    model.eval()
    logger.info('Encoder loaded from %s', path)
    return model, ent2id, rel2id
