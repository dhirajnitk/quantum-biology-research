"""
gnn_pipeline.py  —  Graph Neural Network predictor for FMO quantum transport.

Predicts ENAQT transport efficiency from structural features alone,
bypassing full differential equation execution.

Architecture: 2-layer Graph Convolution Network (GCN)
  Node features: [site energy, connectivity count, dielectric shielding]
  Adjacency: from native Hamiltonian coupling matrix
  Output: predicted transport efficiency

Reference:
  Tang et al. (2025) DeepQT — ML for quantum transport
"""

import numpy as np


class SimplifiedFmoGnnPredictor:
    """Lightweight GCN for FMO transport efficiency prediction."""

    def __init__(self, n_features=3, hidden_dim=8):
        self.W1 = np.random.randn(n_features, hidden_dim) * 0.1
        self.W2 = np.random.randn(hidden_dim, 1) * 0.1

    def relu(self, x):
        return np.maximum(0.0, x)

    def forward_pass(self, adjacency_matrix, node_features):
        """GCN forward: X' = ReLU(A_hat @ X @ W1), then global pool + regression."""
        A_hat = adjacency_matrix + np.eye(adjacency_matrix.shape[0])
        deg = np.sum(A_hat, axis=1)
        D_inv = np.diag(1.0 / np.sqrt(np.clip(deg, 1e-6, None)))
        A_norm = D_inv @ A_hat @ D_inv

        H = self.relu(A_norm @ node_features @ self.W1)
        h_graph = np.mean(H, axis=0)
        return float((h_graph @ self.W2)[0])

    def train_step(self, adj, features, target=None, lr=0.01, target_efficiency=None):
        """Single training step via numerical gradient.

        Supports both positional target and keyword target_efficiency.
        """
        tgt = target if target is not None else target_efficiency
        if tgt is None:
            return 0.0
        pred = self.forward_pass(adj, features)
        loss = (pred - tgt) ** 2

        # Numerical gradient
        delta = 1e-5
        for i in range(self.W2.shape[0]):
            self.W2[i, 0] += delta
            pred_new = self.forward_pass(adj, features)
            grad = ((pred_new - tgt) ** 2 - loss) / delta
            self.W2[i, 0] -= delta
            self.W2[i, 0] -= lr * grad

        return loss

    # API compatibility alias
    train_mock_step = train_step
