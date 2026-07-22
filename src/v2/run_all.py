"""
run_all.py — Orchestrate all v2 computational upgrades.

Runs each of the 5 upgraded analyses in sequence and prints a
summary of results vs the manuscript's current values.
"""

import sys, os, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))


def run_heom():
    print("\n" + "="*60)
    print("UPGRADE 4: HEOM Non-Markovian Dynamics")
    print("="*60)
    from src.v2.heom_dynamics import run_heom, compare_lindblad

    for pid in ['1BL8', '6PV7']:
        print(f"\n--- {pid} ---")
        result = run_heom(pid, t_max_ps=2.0, n_steps=100)
        if result:
            final = result['populations'][:, -1]
            print(f"  Sites: {result['n_sites']}")
            print(f"  Final populations: {[f'{x:.3f}' for x in final]}")


def run_raytrace():
    print("\n" + "="*60)
    print("UPGRADE 5: 3D Ray-Tracing for CCO Path Loss")
    print("="*60)
    from src.v2.ray_tracing_3d import simulate_photon_burst, ensemble_from_raytrace

    result = simulate_photon_burst(n_photons=50000)
    p = result['p_hit']
    print(f"  Photons: {result['n_photons']}")
    print(f"  p_hit (ray-trace) = {p:.2e}")
    print(f"  Paper p = 7.84e-4")
    print(f"  P_success(5000) ray-trace = {ensemble_from_raytrace(p):.4%}")
    print(f"  P_success(5000) paper = {ensemble_from_raytrace(7.84e-4):.4%}")


def run_sdp():
    print("\n" + "="*60)
    print("UPGRADE 3: SDP KCKAS Contextuality Proof")
    print("="*60)
    import numpy as np
    from src.v2.kckas_sdp import kckas_projectors, compute_s_max, compute_fcoh_from_hamiltonian, compute_s_heuristic
    from src.pdb_tools.trp_extractor import fetch_pdb, extract_trp_data, distance_matrix, build_hamiltonian

    try:
        import cvxpy
        P = kckas_projectors()
        S_sdp, _ = compute_s_max(P)
        if S_sdp is not None:
            print(f"  True quantum bound S_max = {S_sdp:.4f} (theoretical: sqrt(5)={np.sqrt(5):.4f})")
            print(f"  Classical bound = 2.0")
            print(f"  S_max > 2: {S_sdp > 2.0}")
    except Exception as e:
        print(f"  SDP not available: {e}")

    # Compare heuristic for targets
    print("\n  Heuristic scores:")
    for pid in ['1BL8', '6PV7', '7TYO']:
        text = fetch_pdb(pid, cache=True)
        if not text:
            continue
        centroids, dipoles = extract_trp_data(text)
        D, keys = distance_matrix(centroids)
        H_cm, _ = build_hamiltonian(D, keys, dipoles, dielectric=2.0)
        F = compute_fcoh_from_hamiltonian(H_cm)
        S = compute_s_heuristic(P, F)
        print(f"    {pid}: S = {S:.4f}")

    import numpy as np


def run_polaritonic():
    print("\n" + "="*60)
    print("UPGRADE 2: Polaritonic Clock Derivation")
    print("="*60)
    from src.v2.polaritonic_clock import derive_effective_hamiltonian, compute_effective_coupling

    result = derive_effective_hamiltonian()
    print(f"  Single-pair J_eff = {result['J_eff_numeric_cm']:.1f} cm^-1")
    print()

    print("  Collective enhancement:")
    for N, g in [(5, 50), (10, 50), (10, 80)]:
        J = compute_effective_coupling(N_sites=N, g_cm=g)
        print(f"    N={N}, g={g} cm^-1: J_eff = {J:.0f} cm^-1")
    print(f"\n  Paper requires: J_max = 400 cm^-1")


def run_transition_dipoles():
    print("\n" + "="*60)
    print("UPGRADE 1: Transition Dipole Analysis (L_a/L_b)")
    print("="*60)
    from src.v2.transition_dipoles import LITERATURE_DIPOLES, compute_full_coupling_matrix

    print("  Literature values for indole:")
    for state, props in LITERATURE_DIPOLES.items():
        print(f"    {state}: {props['wavelength_nm']} nm, "
              f"f = {props['oscillator_strength']}, "
              f"direction = {props['direction']}")

    print("\n  Full coupling matrix (L_a + L_b) for KcsA:")
    result = compute_full_coupling_matrix('1BL8')
    if result is not None:
        J_full, keys = result
        n = len(keys)
        max_aa = np.max(np.abs(J_full[:n, :n]))
        max_bb = np.max(np.abs(J_full[n:, n:]))
        max_ab = np.max(np.abs(J_full[:n, n:]))
        print(f"    Max L_a-L_a: {max_aa:.1f} cm^-1")
        print(f"    Max L_b-L_b: {max_bb:.1f} cm^-1")
        print(f"    Max L_a-L_b cross: {max_ab:.1f} cm^-1")


if __name__ == '__main__':
    t0 = time.time()

    run_transition_dipoles()
    run_polaritonic()
    run_sdp()
    run_heom()
    run_raytrace()

    t1 = time.time()
    print(f"\n{'='*60}")
    print(f"Total time: {t1-t0:.1f} seconds")
    print(f"{'='*60}")
