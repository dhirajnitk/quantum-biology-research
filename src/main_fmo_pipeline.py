"""
main_fmo_pipeline.py  —  Master integration script for the FMO quantum biology pipeline.

Runs all modules sequentially:
  1. Lindblad ENAQT transport
  2. Quantum Darwinism analysis
  3. Holevo capacity + entropy production
  4. GNN transport prediction

Outputs manuscript-ready data for papers P1, P2, P3.
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from src.core.lindblad_solver import FmoLindbladSolver
from src.analysis.quantum_darwinism import FmoDarwinismAnalyzer
from src.analysis.channel_capacity import FmoThermodynamicEngine
from src.ml.gnn_pipeline import SimplifiedFmoGnnPredictor


def execute_complete_fmo_pipeline():
    print("=" * 60)
    print("  SUB-TUBULIN & FMO QUANTUM BIOLOGY BACKEND")
    print("  ENAQT + Darwinism + Capacity + ML Pipeline")
    print("=" * 60)

    # 1. Lindblad Transport
    print("\n[1/4] Lindblad Open Systems Solver...")
    solver = FmoLindbladSolver(disorder_std=50.0)  # 50 cm^-1 disorder for ENAQT
    efficiency, density_history = solver.run_time_evolution(dephasing_rate=2.0)
    print(f"      Transport efficiency: {efficiency*100:.2f}%")

    # ENAQT curve
    print("\n      ENAQT Efficiency vs Dephasing:")
    for g in [0, 10, 50, 100, 175, 300, 500]:
        eff, _ = solver.run_time_evolution(dephasing_rate=g * 0.0188)
        print(f"        gamma={g:>4d} cm-1 -> {eff*100:5.2f}%")

    # QMI at optimal dephasing
    rho_last = density_history[-1]
    qmi = solver.calculate_quantum_mutual_information(rho_last)
    print("\n      QMI Matrix (bits) at t=10ps:")
    for i in range(7):
        row = '  '.join(f"{qmi[i,j]:.4f}" for j in range(7))
        print(f"        {row}")

    # 2. Quantum Darwinism
    print("\n[2/4] Quantum Darwinism Analysis...")
    darwin = FmoDarwinismAnalyzer(density_history)
    darwin.print_report()

    # 3. Thermodynamics
    print("\n[3/4] Thermodynamic & Channel Capacity...")
    H_ps = solver.H_base * 0.0188  # convert cm^-1 to ps^-1
    thermo = FmoThermodynamicEngine(density_history, H_ps)
    holevo = thermo.compute_holevo_capacity(time_index=50)
    entropy_prod = thermo.compute_entropy_production(time_index=50)
    print(f"      Holevo capacity at t=0.5ps:     {holevo:.4f} bits")
    print(f"      Entropy production at t=0.5ps:  {entropy_prod:.4f} kB/ps")

    # 4. GNN Prediction
    print("\n[4/4] Geometric Graph Neural Network...")
    adj = (np.abs(solver.H_base) > 0.1).astype(float)
    np.fill_diagonal(adj, 0)
    features = np.random.randn(7, 3)
    features[:, 0] = np.diag(solver.H_base) / 500.0

    gnn = SimplifiedFmoGnnPredictor()
    loss = gnn.train_mock_step(adj, features, target_efficiency=efficiency)
    predicted = gnn.forward_pass(adj, features)
    print(f"      Target efficiency:    {efficiency:.4f}")
    print(f"      GNN prediction:       {predicted:.4f}")
    print(f"      Training loss:        {loss:.6f}")

    print("\n" + "=" * 60)
    print("  ALL MODULES COMPLETED — DATA READY FOR P1/P2/P3")
    print("=" * 60)


if __name__ == "__main__":
    execute_complete_fmo_pipeline()
