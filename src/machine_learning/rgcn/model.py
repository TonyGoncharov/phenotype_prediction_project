"""
src/machine_learning/rgcn/model.py - R-GCN encoder for heterogeneous gene knowledge graphs.

Implements a two-layer Relational Graph Convolutional Network (R-GCN) that
produces one embedding vector per entity. Used as a gene embedding extractor:
the encoder is trained via link reconstruction on context edges (GO, expression,
PPI, gene-protein), then frozen gene embeddings are passed to a downstream
multi-label classifier for phenotype-class prediction.
"""

from __future__ import annotations

import torch
import torch.nn.functional as F
from torch_geometric.nn import RGCNConv


class RGCNEncoder(torch.nn.Module):
    """
    Two-layer R-GCN that produces one embedding vector per entity.

    The forward pass runs message passing over a heterogeneous edge set
    (bidirectional: forward + inverse relation types). Entity embeddings
    are initialised with Xavier uniform and updated through two RGCNConv
    layers with a ReLU activation between them.

    Args:
        num_entities: Total number of entities in the context graph.
        num_relations: Total number of relation types including inverses
            (i.e. 2 × n_forward_relations for a bidirectional graph).
        hidden_dim: Embedding dimensionality (same for input and output).
    """

    def __init__(self, num_entities: int, num_relations: int, hidden_dim: int) -> None:
        """
        Initialise entity embeddings and two RGCNConv layers.

        Args:
            num_entities: Total number of entities in the context graph.
            num_relations: Total relation types including inverses.
            hidden_dim: Embedding dimensionality used throughout both layers.
        """
        super().__init__()
        self.entity_emb = torch.nn.Embedding(num_entities, hidden_dim)
        self.conv1 = RGCNConv(hidden_dim, hidden_dim, num_relations)
        self.conv2 = RGCNConv(hidden_dim, hidden_dim, num_relations)
        torch.nn.init.xavier_uniform_(self.entity_emb.weight)

    def forward(
        self,
        edge_index: torch.Tensor,
        edge_type: torch.Tensor,
    ) -> torch.Tensor:
        """
        Run two-layer R-GCN message passing over all entities.

        Args:
            edge_index: Edge connectivity tensor, shape (2, num_edges).
            edge_type: Relation type indices, shape (num_edges,).

        Returns:
            Entity embedding matrix, shape (num_entities, hidden_dim).
        """
        x = self.entity_emb.weight
        x = F.relu(self.conv1(x, edge_index, edge_type))
        x = self.conv2(x, edge_index, edge_type)
        return x

    def decode(
        self,
        z: torch.Tensor,
        head: torch.Tensor,
        tail: torch.Tensor,
    ) -> torch.Tensor:
        """
        Dot-product link score used for link-reconstruction training.

        Args:
            z: Entity embedding matrix from forward(), shape (num_entities, hidden_dim).
            head: Head entity indices, shape (batch,).
            tail: Tail entity indices, shape (batch,).

        Returns:
            Scores tensor, shape (batch,).
        """
        return (z[head] * z[tail]).sum(dim=-1)
