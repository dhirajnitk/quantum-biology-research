"""
run_tests.py  —  Comprehensive verification of all modules.

Tests:
  1. PDB extraction + dipole orientation
  2. Hamiltonian builder with kappa factor
  3. Biophoton relay + Z-channel
  4. KCKAS contextuality (static + clock-driven)
  5. Lindblad solver (ENAQT + QMI)
  6. Quantum Darwinism (redundancy, pointer states, SBS)
  7. Thermodynamics (Holevo, entropy production)
  8. GNN transport predictor (train + evaluate)
"""

import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

n_pass = 0
n_fail = 0

def test(name, condition, detail=""):
    global n_pass, n_fail
    if condition:
        n_pass += 1
        print(f"  [PASS] {name}")
    else:
        n_fail += 1
        print(f"  [FAIL] {name} {detail}")


def run_all():
    global n_pass, n_fail
    print("=" * 68)
    print("  COMPREHENSIVE TEST SUITE")
    print("  Sub-Tubular Quantum Information Processing")
    print("=" * 68)

    # ── 1. PDB extraction + Dipole orientation ────────────────
    print("\n[1] PDB TOOLS + DIPOLE ORIENTATION")
    from src.pdb_tools.trp_extractor import (
        fetch_pdb, extract_trp_data, extract_trp_coordinates,
        distance_matrix, dipole_orientation_factor, build_hamiltonian
    )
    text = fetch_pdb("1BL8")
    test("fetch_pdb(1BL8)", text is not None)

    centroids, dipoles = extract_trp_data(text)
    test("extract_trp_data centroids", len(centroids) == 5)
    test("extract_trp_data dipoles", len(dipoles) == 5)
    for seq, d in dipoles.items():
        test(f"dipole unit vector Trp-{seq}", abs(np.linalg.norm(d) - 1.0) < 1e-6,
             f"norm={np.linalg.norm(d):.4f}")

    D, keys = distance_matrix(centroids)
    test("distance_matrix shape", D.shape == (5, 5))

    # Kappa factor
    mu_i = np.array([1.0, 0.0, 0.0])
    mu_j = np.array([0.0, 1.0, 0.0])
    R = np.array([1.0, 0.0, 0.0])
    k = dipole_orientation_factor(mu_i, mu_j, R)
    test("dipole orientation (orthogonal=0)", abs(k) < 1e-10)

    mu_j = np.array([1.0, 0.0, 0.0])
    k = dipole_orientation_factor(mu_i, mu_j, R)
    test("dipole orientation (collinear=-2)", abs(k - (-2.0)) < 1e-10)

    H, hk = build_hamiltonian(D, keys, dipoles, centroids=centroids)
    test("Hamiltonian builder with dipoles", H.shape == (5, 5))

    # ── 2. Core Hamiltonian ───────────────────────────────────
    print("\n[2] CORE HAMILTONIAN")
    from src.core.hamiltonian import (
        build_pdb_hamiltonian, build_clock_hamiltonian,
        build_multiscale_hamiltonian, coherence_factor, kckas_s,
        calculate_orientation_factor
    )
    coords = np.array([centroids[k] for k in sorted(centroids)])
    H_pdb, D_out = build_pdb_hamiltonian(coords)
    test("build_pdb_hamiltonian", H_pdb.shape == (5, 5))

    k2 = calculate_orientation_factor(mu_i, mu_j, R)
    test("calculate_orientation_factor", abs(k2 - (-2.0)) < 1e-10)

    H_total, _ = build_multiscale_hamiltonian(coords, clock_amplitude=1.0)
    test("build_multiscale_hamiltonian", H_total.shape == (5, 5))

    F = coherence_factor(H_total)
    test("coherence_factor range", 0 <= F <= 1, f"F={F:.4f}")

    S = kckas_s(H_total)
    test("kckas_s range", 0 <= S <= 2.2361, f"S={S:.4f}")

    # ── 3. Biophoton relay + Z-channel ────────────────────────
    print("\n[3] BIOPHOTON RELAY")
    from src.core.biophoton_relay import BiophotonRelay, z_channel_capacity, TRP_QY
    test("TRP_QY", abs(TRP_QY - 0.21) < 0.01)

    relay = BiophotonRelay(2.5)
    for target in ["trp", "cco"]:
        p_ens, p_hit, n_arr = relay.spatial_ensemble_success(1000, target_type=target)
        test(f"spatial ensemble {target} N=1000", 0 < p_ens < 1, f"p={p_ens:.4f}")

    p_ens, p_hit, _ = relay.spatial_ensemble_success(5000, target_type="cco")
    test("CCO N=5000 > 90%", p_ens > 0.9, f"p={p_ens:.4f}")

    p_ens, p_hit, _ = relay.spatial_ensemble_success(10000, target_type="cco")
    test("CCO N=10000 > 99%", p_ens > 0.99, f"p={p_ens:.4f}")

    c = z_channel_capacity(0.95)
    test("z_channel_capacity(0.95)", 0 < c <= 1, f"C={c:.4f}")

    # ── 4. KCKAS contextuality ────────────────────────────────
    print("\n[4] KCKAS CONTEXTUALITY")
    from src.analysis.pdb_contextuality import compute_kckas_from_pdb, KCKAS_CLASSICAL_BOUND
    from src.analysis.quantum_hamiltonian_engine import HamiltonianEngine

    r = compute_kckas_from_pdb("1BL8")
    test("KCKAS 1BL8", r is not None, f"S={r['kckas_sum']:.4f}")
    if r:
        test("KCKAS 1BL8 S < 2", r["kckas_sum"] < 2.0, f"S={r['kckas_sum']:.4f}")

    eng = HamiltonianEngine("1BL8")
    S0 = eng.compute_kckas_s(0.0)
    S1 = eng.compute_kckas_s(1.0)
    test("Hamiltonian S(0) matches KCKAS", abs(S0 - r["kckas_sum"]) < 0.01, f"S0={S0:.4f}")
    test("Hamiltonian S(1) > 2", S1 > 2.0, f"S1={S1:.4f}")

    # ── 5. Lindblad solver ────────────────────────────────────
    print("\n[5] LINDBLAD SOLVER")
    from src.core.lindblad_solver import FmoLindbladSolver

    solver = FmoLindbladSolver(disorder_std=50.0)
    eff, hist = solver.run_time_evolution(dephasing_rate=175 * 0.0188)
    test("Lindblad evolution", len(hist) > 0, f"steps={len(hist)}, eff={eff*100:.1f}%")

    qmi = solver.compute_qmi_matrix(hist[-1])
    test("QMI matrix shape", qmi.shape == (7, 7))

    # ENAQT curve
    gammas, effs = solver.enaqt_curve()
    test("ENAQT curve", len(effs) == len(gammas))
    # Check there's a non-monotonic ENAQT peak
    is_enaqt = any(effs[i] > effs[i-1] and effs[i] > effs[i+1] for i in range(1, len(effs)-1))
    test("ENAQT non-monotonic (peak)", is_enaqt)

    # ── 6. Thermodynamics ─────────────────────────────────────
    print("\n[6] THERMODYNAMICS")
    from src.analysis.channel_capacity import FmoThermodynamicEngine

    thermo = FmoThermodynamicEngine(hist, solver.H)
    S_t = thermo.compute_entropy_trajectory()
    test("entropy trajectory", len(S_t) == len(hist))
    test("entropy > 0", np.max(S_t) > 0)

    chi = thermo.compute_holevo_capacity(time_index=50)
    test("Holevo capacity > 0", chi > 0, f"chi={chi:.4f}")

    dS = thermo.compute_entropy_production(time_index=50)
    test("entropy production > 0", dS >= 0, f"dS={dS:.6f}")

    # ── 7. Quantum Darwinism ──────────────────────────────────
    print("\n[7] QUANTUM DARWINISM")
    from src.analysis.quantum_darwinism import QuantumDarwinism

    qd = QuantumDarwinism(hist, system_site=0)
    profile = qd.compute_redundancy_curve(time_index=50)
    test("redundancy curve", len(profile) == 6)

    mi = qd.mi_between_sites(0, 1, time_index=50)
    test("pairwise QMI > 0", mi > 0, f"I(1:2)={mi:.4f}")

    ps, vals = qd.find_pointer_states(time_index=100)
    test("pointer states", len(vals) > 0, f"n_evals={len(vals)}")

    sbs = qd.compute_sbs_fidelity(time_index=100)
    test("SBS fidelity", 0 <= sbs <= 1, f"SBS={sbs:.4f}")

    # ── 8. GNN ────────────────────────────────────────────────
    print("\n[8] GNN TRANSPORT PREDICTOR")
    from src.ml.gnn_transport_predictor import StructuralQuantumGnn, QuantumTrainingDataFactory

    gnn = StructuralQuantumGnn(input_dim=4, hidden_dim=12)
    factory = QuantumTrainingDataFactory()

    A, X, targets = factory.construct_synthetic_protein_node()
    pred = gnn.forward(A, X)
    test("GNN forward", len(pred) == 2)

    loss = gnn.train_step(A, X, targets)
    test("GNN train_step", loss >= 0)

    # Quick training
    losses = []
    for _ in range(200):
        A2, X2, t2 = factory.construct_synthetic_protein_node()
        losses.append(gnn.train_step(A2, X2, t2, lr=0.001))
    test("GNN training converges", sum(losses[-50:]) / 50 < 10.0,
         f"final_loss={sum(losses[-50:])/50:.4f}")

    # ── Summary ───────────────────────────────────────────────
    print(f"\n{'=' * 68}")
    print(f"  RESULTS: {n_pass} passed, {n_fail} failed out of {n_pass + n_fail} tests")
    if n_fail == 0:
        print("  ALL TESTS PASSED")
    else:
        print(f"  {n_fail} TEST(S) FAILED")
    print(f"{'=' * 68}")
    return n_fail == 0


if __name__ == "__main__":
    success = run_all()
    sys.exit(0 if success else 1)
