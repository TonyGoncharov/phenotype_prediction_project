from __future__ import annotations

import torch
import torch.nn.functional as F
from torch_geometric.nn import RGCNConv


class RGCNEncoder(torch.nn.Module):
    """Two-layer R-GCN that produces one embedding vector per entity."""

    def __init__(self, num_entities: int, num_relations: int, hidden_dim: int):
        super().__init__()
        self.entity_emb = torch.nn.Embedding(num_entities, hidden_dim)
        self.conv1 = RGCNConv(hidden_dim, hidden_dim, num_relations)
        self.conv2 = RGCNConv(hidden_dim, hidden_dim, num_relations)
        torch.nn.init.xavier_uniform_(self.entity_emb.weight)

    def forward(self, edge_index: torch.Tensor, edge_type: torch.Tensor) -> torch.Tensor:
        x = self.entity_emb.weight
        x = F.relu(self.conv1(x, edge_index, edge_type))
        x = self.conv2(x, edge_index, edge_type)
        return x

    def decode(self, z: torch.Tensor, head: torch.Tensor, tail: torch.Tensor) -> torch.Tensor:
        """Dot-product score (used for link-reconstruction training)."""
        return (z[head] * z[tail]).sum(dim=-1)
