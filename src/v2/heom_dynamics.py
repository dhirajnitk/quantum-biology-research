"""
heom_dynamics.py — Non-Markovian HEOM trajectories for Trp networks.

QuTiP 5 includes a HEOM solver that captures non-Markovian memory effects
missed by the Lindblad master equation. We run full HEOM trajectories on
the top Trp networks (KcsA, nAChR, Nav1.4) and compare population dynamics
with the Markovian Lindblad result.
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

# Conditional import — HEOM solver requires QuTiP 5
try:
    import qutip as qt
    from qutip.solver.heom import HEOMSolver
    from qutip.solver.heom.bofin_baths import DrudeLorentzBath
    if qt.__version__.startswith('5'):
        HAS_QUTIP = True
    else:
        HAS_QUTIP = False
except Exception:
    HAS_QUTIP = False

from src.pdb_tools.trp_extractor import fetch_pdb, extract_trp_data, distance_matrix

# Constants
HBAR = 1.054571817e-34       # J·s
KB = 1.380649e-23             # J/K
CM1_TO_J = 1.9863e-23         # J per cm^-1
CM1_TO_RADS = 2 * np.pi * 2.99792458e10  # rad/s per cm^-1


def build_trp_hamiltonian(pdb_id, eps=2.0, disorder_std=30.0):
    """Build tight-binding Hamiltonian from PDB Trp coordinates."""
    text = fetch_pdb(pdb_id, cache=True)
    if not text:
        return None, None
    centroids, _ = extract_trp_data(text)
    keys = list(centroids.keys())
    if len(keys) < 3:
        return None, None

    n = len(keys)
    D, _ = distance_matrix(centroids)

    J0 = -80.0
    R0 = 10.0
    H_cm = np.zeros((n, n))
    for i in range(n):
        H_cm[i, i] = float(i * 80) + np.random.normal(0, disorder_std)
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] > 0.1:
                J_ij = J0 * (R0 / D[i, j]) ** 3 / np.sqrt(eps)
                H_cm[i, j] = H_cm[j, i] = J_ij

    H = qt.Qobj(H_cm * CM1_TO_RADS)
    return H, keys


def run_heom(pdb_id, temperature=310.0, t_max_ps=5.0, n_steps=500):
    """Run HEOM dynamics on a Trp network and return trajectories."""
    if not HAS_QUTIP:
        print("HEOM solver requires QuTiP 5.")
        return None

    H, keys = build_trp_hamiltonian(pdb_id)
    if H is None:
        print(f"Could not build Hamiltonian for {pdb_id}")
        return None

    n = H.shape[0]

    # Initial state: excitation on site 0
    rho0 = qt.ket2dm(qt.basis(n, 0))

    # Bath parameters: Drude-Lorentz spectral density
    kB = 1.380649e-23
    hbar = 1.054571817e-34
    lam = 35.0 * 1.9863e-23      # reorganisation energy (J)
    gamma_c = 53.0 * 2 * np.pi * 2.99792458e10  # cutoff frequency (rad/s)
    kT = kB * temperature

    # System-bath coupling operators: dephasing on each site
    Q = [qt.basis(n, i) * qt.basis(n, i).dag() for i in range(n)]

    # Create Drude-Lorentz baths
    baths = []
    for q in Q:
        bath = DrudeLorentzBath(q, lam, gamma_c, kT, Nk=3)
        baths.append(bath)

    # HEOM solver
    tlist = np.linspace(0, t_max_ps * 1e-12, n_steps)

    try:
        solver = HEOMSolver(H, baths, max_depth=5)

        # Expectation operators: site populations
        e_ops = [qt.basis(n, i) * qt.basis(n, i).dag() for i in range(n)]

        result = solver.run(rho0, tlist, e_ops=e_ops)
        pops = np.array(result.expect)
        times = result.times

        return {
            'pdb_id': pdb_id,
            'times': times,
            'populations': pops,
            'n_sites': n,
            'keys': keys,
        }
    except Exception as e:
        print(f"HEOM error for {pdb_id}: {e}")
        return None


def compare_lindblad(pdb_id, temperature=310.0, t_max_ps=5.0, n_steps=500):
    """Compare HEOM vs Lindblad population dynamics."""
    heom_result = run_heom(pdb_id, temperature, t_max_ps, n_steps)
    if heom_result is None:
        return None

    H, _ = build_trp_hamiltonian(pdb_id)
    n = H.shape[0]
    rho0 = qt.ket2dm(qt.basis(n, 0))
    tlist = np.linspace(0, t_max_ps * 1e-12, n_steps)

    # Lindblad dephasing rate (same as HEOM bath)
    kT = KB * temperature
    lam = 35.0 * CM1_TO_J
    gamma_c = 53.0 * CM1_TO_RADS
    gamma_phi = (2 * np.pi * kT / HBAR) * (lam / gamma_c)

    c_ops = [np.sqrt(gamma_phi) * qt.basis(n, i) * qt.basis(n, i).dag() for i in range(n)]
    e_ops = [qt.basis(n, i) * qt.basis(n, i).dag() for i in range(n)]

    result = qt.mesolve(H, rho0, tlist, c_ops, e_ops)
    lindblad_pops = np.array(result.expect)

    return {
        'heom': heom_result,
        'lindblad_times': tlist,
        'lindblad_populations': lindblad_pops,
    }


if __name__ == '__main__':
    for pid in ['1BL8', '6PV7', '7TYO']:
        print(f"\n=== HEOM: {pid} ===")
        result = run_heom(pid, t_max_ps=2.0, n_steps=200)
        if result:
            final_pops = result['populations'][:, -1]
            print(f"  Sites: {result['n_sites']}")
            print(f"  Final populations: {np.round(final_pops, 4)}")
            print(f"  Coherence fraction: {1 - np.sum(final_pops**2):.4f}")
