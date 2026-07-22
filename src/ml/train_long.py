"""
train_long.py — Extended ML training pipeline (runs hours, saves model).
Generates 20k+ samples, tunes hyperparameters, saves best model.
"""
import numpy as np, sys, os, pickle, time, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from src.ml.fast_transport_predictor import generate_dataset

print("=" * 60)
print("LONG-RUNNING ML TRAINING PIPELINE")
print("=" * 60)

# Phase 1: Generate 20000 samples
print("\n[Phase 1] Generating 20000 structures...")
t0 = time.time()
X, y = generate_dataset(20000, seed=42)
print(f"  Generated {len(X)} samples in {time.time()-t0:.1f}s")
print(f"  Features: {X.shape[1]}, Targets: {y.shape[1]}")
print(f"  ENAQT range: [{y[:,0].min():.4f}, {y[:,0].max():.4f}]")
print(f"  PropDelay range: [{y[:,2].min():.4f}, {y[:,2].max():.4f}] ps")

# Split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.15, random_state=42)

# Phase 2: RandomForest with hyperparameter search
print("\n[Phase 2] RandomForest hyperparameter search...")
param_grid = {
    'n_estimators': [100, 200, 500, 1000],
    'max_depth': [5, 10, 15, 20, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
}

target_names = ['ENAQT', 'Holevo', 'PropDelay']
best_models = {}
results = {}

for i in range(3):
    print(f"\n  Training {target_names[i]}...")
    t0 = time.time()
    
    rf = RandomForestRegressor(random_state=42, n_jobs=-1, verbose=0)
    search = RandomizedSearchCV(rf, param_grid, n_iter=20, cv=3, 
                                 scoring='r2', n_jobs=-1, random_state=42)
    search.fit(X_train, y_train[:, i])
    
    best_model = search.best_estimator_
    y_pred = best_model.predict(X_test)
    
    ss_res = np.sum((y_test[:, i] - y_pred)**2)
    ss_tot = np.sum((y_test[:, i] - y_test[:, i].mean())**2)
    r2 = 1 - ss_res / max(ss_tot, 1e-10)
    mae = np.mean(np.abs(y_test[:, i] - y_pred))
    
    print(f"  Best params: {search.best_params_}")
    print(f"  R² = {r2:.4f}, MAE = {mae:.4f}")
    print(f"  Time: {time.time()-t0:.1f}s")
    
    best_models[target_names[i]] = best_model
    results[target_names[i]] = {'r2': r2, 'mae': mae, 'best_params': search.best_params_}

# Phase 3: PyTorch MLP training (if torch available)
try:
    import torch
    import torch.nn as nn
    print("\n[Phase 3] PyTorch MLP training...")
    
    scaler_X = StandardScaler()
    X_train_s = scaler_X.fit_transform(X_train)
    X_test_s = scaler_X.transform(X_test)
    
    class MLP(nn.Module):
        def __init__(self, input_dim=9):
            super().__init__()
            self.net = nn.Sequential(
                nn.Linear(input_dim, 128), nn.ReLU(), nn.Dropout(0.1),
                nn.Linear(128, 64), nn.ReLU(), nn.Dropout(0.1),
                nn.Linear(64, 32), nn.ReLU(),
                nn.Linear(32, 3), nn.Sigmoid())
        def forward(self, x):
            return self.net(x)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = MLP().to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.001, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=2000, eta_min=1e-6)
    
    X_train_t = torch.tensor(X_train_s, dtype=torch.float32).to(device)
    y_train_t = torch.tensor(y_train, dtype=torch.float32).to(device)
    X_test_t = torch.tensor(X_test_s, dtype=torch.float32).to(device)
    y_test_t = torch.tensor(y_test, dtype=torch.float32).to(device)
    
    best_test_loss = float('inf')
    for epoch in range(2000):
        model.train()
        optimizer.zero_grad()
        pred = model(X_train_t)
        loss = nn.MSELoss()(pred, y_train_t)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        optimizer.step()
        scheduler.step()
        
        if epoch % 200 == 0:
            model.eval()
            with torch.no_grad():
                p_test = model(X_test_t)
                test_loss = nn.MSELoss()(p_test, y_test_t).item()
                if test_loss < best_test_loss:
                    best_test_loss = test_loss
                    torch.save(model.state_dict(), 'src/ml/best_mlp.pt')
                r2_e = 1 - torch.sum((y_test_t[:,0]-p_test[:,0])**2) / max(torch.sum((y_test_t[:,0]-y_test_t[:,0].mean())**2), 1)
                print(f"  Epoch {epoch:4d}: train={loss.item():.6f} test={test_loss:.6f} R²_e={r2_e:.4f}")
    
    results['MLP_ENAQT_R2'] = r2_e.item()
    results['MLP_test_loss'] = best_test_loss
    
except Exception as e:
    print(f"  PyTorch skipped: {e}")

# Save all results
with open('src/ml/training_results.pkl', 'wb') as f:
    pickle.dump({'results': results, 'models': best_models}, f)

print("\n" + "=" * 60)
print("TRAINING COMPLETE")
print(json.dumps(results, indent=2, default=str))
print("=" * 60)
