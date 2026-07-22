"""
fast_transport_predictor.py — Fast quantum transport prediction from graph features.

Strategy:
1. Generate diverse chromophore geometries (N=4-12, varied positions/dipoles/dephasing)
2. Extract graph-level features that capture the relevant physics
3. Train PyTorch MLP to predict ENAQT efficiency and Holevo capacity
4. Train/validate split to prevent overfitting
"""

import torch
import torch.nn as nn
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


def compute_coupling(mu_i, mu_j, R_ij, C=80.0):
    r = np.linalg.norm(R_ij)
    if r < 0.1:
        return 0.0
    R_hat = R_ij / r
    kappa = np.dot(mu_i, mu_j) - 3 * np.dot(mu_i, R_hat) * np.dot(mu_j, R_hat)
    return C * kappa / (r ** 3)


def generate_dataset(n_samples=2000, seed=42):
    """Generate dataset with graph-level features and physics-based targets."""
    np.random.seed(seed)
    features = []
    targets = []

    for _ in range(n_samples):
        n = np.random.randint(4, 13)
        target_site = n - 1

        # Random 3D positions
        radius = np.random.uniform(10, 35)
        pos = np.random.normal(0, radius / 3, (n, 3))
        pos -= pos.mean(axis=0)

        dip = np.random.normal(0, 1, (n, 3))
        dip = dip / (np.linalg.norm(dip, axis=1, keepdims=True) + 1e-10)

        deph = np.random.uniform(1, 300)

        # Build coupling matrix
        couplings = []
        for i in range(n):
            for j in range(i + 1, n):
                c = compute_coupling(dip[i], dip[j], pos[i] - pos[j])
                couplings.append(c)
        couplings = np.array(couplings)

        # Graph-level features (9-dim)
        mean_coupling = np.mean(np.abs(couplings)) if len(couplings) > 0 else 0
        max_coupling = np.max(np.abs(couplings)) if len(couplings) > 0 else 0
        coupling_spread = np.std(couplings) if len(couplings) > 0 else 0
        n_sites = n
        dip_corr = np.mean([np.abs(np.dot(dip[i], dip[j])) for i in range(n) for j in range(i + 1, n)])
        sink_dist = np.min(np.linalg.norm(pos - pos[target_site], axis=1))
        mean_dist = np.mean([np.linalg.norm(pos[i] - pos[j]) for i in range(n) for j in range(i + 1, n)])
        deph_norm = deph / 200.0
        coupling_deph_ratio = max_coupling / max(deph / 200.0, 0.01)

        feat = np.array([
            mean_coupling, max_coupling, coupling_spread, n_sites / 12.0,
            dip_corr, sink_dist / 35.0, mean_dist / 35.0, deph_norm,
            coupling_deph_ratio,
        ])

        # Physics-based targets: ENAQT efficiency, Holevo capacity, and propagation delay
        J_eff = mean_coupling
        eta = 0.8 * (1 - np.exp(-J_eff / max(deph, 1))) * np.exp(-deph / 120)
        eta = max(0.02, min(0.99, eta))
        holevo = eta * (1 + 0.3 * J_eff / (J_eff + deph + 1))
        holevo = max(0.02, min(2.0, holevo))
        
        # Propagation delay: time to reach 63% of max efficiency (1/e timescale)
        # Proportional to inverse mean coupling, modulated by dephasing
        tau_prop = (1.0 / max(J_eff, 0.01)) * (1 + deph / 100.0)
        tau_prop = min(max(tau_prop, 0.01), 10.0)  # clamp to [0.01, 10] ps

        features.append(feat)
        targets.append([eta, holevo, tau_prop])

    return np.array(features), np.array(targets)


class QuantumMLP(nn.Module):
    """Simple MLP for quantum transport prediction from graph features."""

    def __init__(self, input_dim=9, hidden_dim=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 2),
            nn.Sigmoid(),
        )

    def forward(self, x):
        return self.net(x)


def train_model(model, X_train, y_train, X_val, y_val, epochs=1000, lr=0.001):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs, eta_min=1e-6)

    best_val_loss = float('inf')
    patience = 100
    patience_counter = 0

    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        pred = model(X_train)
        loss = nn.MSELoss()(pred, y_train)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        optimizer.step()
        scheduler.step()

        if epoch % 100 == 0:
            model.eval()
            with torch.no_grad():
                pred_val = model(X_val)
                val_loss = nn.MSELoss()(pred_val, y_val).item()
                train_loss = loss.item()
                mae_e = torch.abs(pred_val[:, 0] - y_val[:, 0]).mean().item()
                ss_res = torch.sum((y_val[:, 0] - pred_val[:, 0]) ** 2).item()
                ss_tot = torch.sum((y_val[:, 0] - y_val[:, 0].mean()) ** 2).item()
                r2 = 1 - ss_res / max(ss_tot, 1e-10)

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_counter = 0
            else:
                patience_counter += 100
                if patience_counter >= patience:
                    print(f"  Early stopping at epoch {epoch}")
                    break

    return best_val_loss


if __name__ == "__main__":
    print("=" * 60)
    print("  FAST QUANTUM TRANSPORT PREDICTOR")
    print(f"  Device: {device}")
    print("=" * 60)

    # Generate data
    print("\nGenerating 2000 structures...")
    X, y = generate_dataset(2000)
    print(f"  Feature dim: {X.shape[1]}")
    print(f"  ENAQT range: [{y[:,0].min():.3f}, {y[:,0].max():.3f}]")
    print(f"  Holevo range: [{y[:,1].min():.3f}, {y[:,1].max():.3f}]")

    # Train/val/test split (60/20/20)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.25, random_state=42)

    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_val_s = scaler.transform(X_val)
    X_test_s = scaler.transform(X_test)

    X_train_t = torch.tensor(X_train_s, dtype=torch.float32).to(device)
    y_train_t = torch.tensor(y_train, dtype=torch.float32).to(device)
    X_val_t = torch.tensor(X_val_s, dtype=torch.float32).to(device)
    y_val_t = torch.tensor(y_val, dtype=torch.float32).to(device)

    print(f"\n  Training: {len(X_train)}, Validation: {len(X_val)}, Test: {len(X_test)}")

    # Train
    model = QuantumMLP(input_dim=9).to(device)
    print("\nTraining...")
    train_model(model, X_train_t, y_train_t, X_val_t, y_val_t)

    # Final test evaluation
    model.eval()
    X_test_t = torch.tensor(X_test_s, dtype=torch.float32).to(device)
    y_test_t = torch.tensor(y_test, dtype=torch.float32).to(device)

    with torch.no_grad():
        pred_test = model(X_test_t)
        mae_e = torch.abs(pred_test[:, 0] - y_test_t[:, 0]).mean().item()
        mae_h = torch.abs(pred_test[:, 1] - y_test_t[:, 1]).mean().item()
        ss_res_e = torch.sum((y_test_t[:, 0] - pred_test[:, 0]) ** 2).item()
        ss_tot_e = torch.sum((y_test_t[:, 0] - y_test_t[:, 0].mean()) ** 2).item()
        r2_e = 1 - ss_res_e / max(ss_tot_e, 1e-10)
        ss_res_h = torch.sum((y_test_t[:, 1] - pred_test[:, 1]) ** 2).item()
        ss_tot_h = torch.sum((y_test_t[:, 1] - y_test_t[:, 1].mean()) ** 2).item()
        r2_h = 1 - ss_res_h / max(ss_tot_h, 1e-10)

        print(f"\n  FINAL TEST RESULTS ({len(X_test)} samples):")
        print(f"    ENAQT MAE:  {mae_e:.4f}  ({mae_e/max(y_test[:,0].mean(),0.01)*100:.1f}%)")
        print(f"    Holevo MAE: {mae_h:.4f}")
        print(f"    ENAQT R²:   {r2_e:.4f}")
        print(f"    Holevo R²:  {r2_h:.4f}")
        print(f"    Pred range: [{pred_test[:,0].min():.3f}, {pred_test[:,0].max():.3f}]")
        print(f"    True range: [{y_test[:,0].min():.3f}, {y_test[:,0].max():.3f}]")
        print("=" * 60)
