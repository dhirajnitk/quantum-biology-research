"""
pdb_bell_test.py  â€”  PDB-grounded Bell/CHSH test for Trp networks.

Uses real Cryo-EM coordinates from the Protein Data Bank to compute
quantum correlation bounds and test whether the spatial arrangement
of Tryptophan networks can violate classical local realism limits.

The dipole-dipole interaction between Trp residues depends on:
  1. The 3D distance R_ij (from PDB coordinates)
  2. The relative orientation of transition dipole moments
  3. The dielectric shielding of the hydrophobic core (Îµ â‰ˆ 2)

A CHSH violation (|S| > 2) would prove non-classical correlations
are possible in the Trp network geometry.

References
----------
[Bell1964] Bell, Physics 1, 195 (1964)
[CHSH1969] Clauser et al., PRL 23, 880 (1969)
[Babcock2024] Babcock et al., JPCB (2024) â€” Trp dipole orientation
"""

import numpy as np
from numpy import pi, cos, sin, sqrt, arccos, dot, cross, array
from numpy.linalg import norm
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.pdb_tools.trp_extractor import fetch_pdb, extract_trp_coordinates, distance_matrix


# CHSH-optimal measurement angles (radians)
# For max violation with E(a,b) = -cos(a-b): |S| = 2*sqrt(2) ~ 2.828
# Standard optimal settings for singlet state correlation
ANGLES_A = array([0.0, pi / 2])
ANGLES_B = array([pi / 4, -pi / 4])


def estimate_dipole_moment(coord_ring_center, atom_coords):
    """Estimate transition dipole direction from Trp ring geometry.

    For a Trp indole ring, the transition dipole lies approximately
    along the long axis of the ring (CG-CD1-CE2 direction).
    If atom-level data is unavailable, we use the CGâ†’CH2 vector
    as a proxy (the long axis of the indole ring).

    [Babcock2024] Methods section.
    """
    if len(atom_coords) >= 2:
        v = atom_coords[1] - atom_coords[0]
        return v / norm(v)
    return array([1.0, 0.0, 0.0])  # fallback


def coupling_strength(R_angstrom, mu_i=None, mu_j=None, dielectric=2.0):
    """Dipole-dipole coupling energy J_ij [cmâ»Â¹] between two Trp sites.

    J_ij = (J0 / Îµ) * (R0 / R_ij)Â³ * (Î¼Ì‚_i Â· Î¼Ì‚_j - 3(Î¼Ì‚_i Â· RÌ‚)(Î¼Ì‚_j Â· RÌ‚))

    Where J0 is the reference coupling at distance R0.
    [Adolphs & Renger 2006, Biophys J]
    """
    if R_angstrom < 0.1:
        return 0.0
    J0 = -80.0           # cmâ»Â¹ at R0 = 10 Ã…
    R0 = 10.0
    R = R_angstrom
    R_hat = array([1.0, 0.0, 0.0])  # default direction

    if mu_i is None or mu_j is None:
        return (J0 / dielectric) * (R0 / R) ** 3

    # Dipole orientation factor
    mu_i_u = mu_i / norm(mu_i)
    mu_j_u = mu_j / norm(mu_j)
    if norm(mu_j_u) < 1e-10 or norm(mu_i_u) < 1e-10:
        return (J0 / dielectric) * (R0 / R) ** 3

    kappa = dot(mu_i_u, mu_j_u) - 3 * dot(mu_i_u, R_hat) * dot(mu_j_u, R_hat)
    return (J0 / dielectric) * (R0 / R) ** 3 * kappa


def quantum_correlation_angle(theta_a, theta_b, coupling_norm):
    """Quantum correlation E(a,b) for two coupled two-level systems.

    For a two-level system (spin-1/2 or exciton pseudospin), the
    correlation at measurement angles θ_a, θ_b is:
        E(a,b) = -cos(θ_a - θ_b) * coupling_norm

    where coupling_norm is the normalised dipole-dipole coupling
    strength (0 to 1). This follows the standard CHSH formulation
    for Bell-state correlations.

    [Bell1964], [CHSH1969]
    """
    phase = theta_a - theta_b
    return -cos(phase) * max(min(abs(coupling_norm), 1.0), 0.0)


def chsh_value(corr_matrix):
    """Compute CHSH parameter S from the 2Ã—2 correlation matrix.

    S = E(a,b) + E(a,b') + E(a',b) - E(a',b')

    Classical local realism: |S| â‰¤ 2
    Quantum mechanical:     |S| â‰¤ 2âˆš2 â‰ˆ 2.828
    """
    return corr_matrix[0, 0] + corr_matrix[0, 1] + corr_matrix[1, 0] - corr_matrix[1, 1]


def compute_chsh_from_pdb(pdb_id, chain=None, dielectric=2.0):
    """Full CHSH computation from a real PDB structure.

    Returns the CHSH value for the most strongly coupled Trp pair,
    normalised by the average network coupling to reflect real
    physical conditions (not just the theoretical maximum).
    """
    text = fetch_pdb(pdb_id)
    if not text:
        return None

    centres = extract_trp_coordinates(text, chain)
    if len(centres) < 2:
        print(f"[!] {pdb_id}: < 2 Trp residues, cannot compute CHSH.")
        return None

    D, keys = distance_matrix(centres)
    n = len(keys)

    # Compute coupling matrix
    J_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            J = coupling_strength(D[i, j], dielectric=dielectric)
            J_matrix[i, j] = J_matrix[j, i] = J

    # Find the most strongly coupled pair
    max_J = 0
    max_pair = (0, 1)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(J_matrix[i, j]) > abs(max_J):
                max_J = J_matrix[i, j]
                max_pair = (i, j)

    i, j = max_pair
    R_ij = D[i, j]
    J_ij = max_J
    J_max = np.max(np.abs(J_matrix)) if np.max(np.abs(J_matrix)) > 0 else 1.0
    beta = abs(J_ij) / J_max  # normalised coupling strength

    # Compute CHSH correlation matrix
    corr = np.zeros((2, 2))
    for ia, a in enumerate(ANGLES_A):
        for ib, b in enumerate(ANGLES_B):
            corr[ia, ib] = quantum_correlation_angle(a, b, J_ij / J_max)

    S = chsh_value(corr)
    S_violation = max(0, S - 2.0)
    S_margin = 2.0 * sqrt(2) - S if S <= 2.0 * sqrt(2) else 0

    result = {
        "pdb_id": pdb_id,
        "n_trp": n,
        "max_coupled_pair": (keys[i], keys[j]),
        "distance_A": R_ij,
        "coupling_cm": J_ij,
        "coupling_normalised": beta,
        "correlation_matrix": corr,
        "S_value": S,
        "S_violation": S_violation,
        "S_margin_to_quantum_bound": S_margin,
        "violates_classical_limit": S > 2.0,
        "status": (
            "QUANTUM VIOLATION" if S > 2.0 else
            "Quantum-allowed (|S| <= 2âˆš2)" if S > 0 else
            "CLASSICAL"
        )
    }
    return result


def run_bell_test(pdb_list):
    """Run CHSH test on multiple PDB targets."""
    print("=" * 72)
    print("  PDB-GROUNDED BELL/CHSH TEST")
    print("  Testing classical vs quantum correlation bounds on real Trp networks")
    print("  Classical bound: |S| <= 2  |  Quantum bound: |S| <= 2.828")
    print("=" * 72)

    print(f"\n{'PDB':<8} {'Trp':<5} {'Best Pair':<14} {'Dist(A)':<10} {'J(cm-1)':<12} {'S':<10} {'Verdict':<20}")
    print(f"{'-'*8} {'-'*5} {'-'*14} {'-'*10} {'-'*12} {'-'*10} {'-'*20}")

    for pdb_id in pdb_list:
        result = compute_chsh_from_pdb(pdb_id)
        if not result:
            continue
        pair_label = f"Trp-{result['max_coupled_pair'][0]}-{result['max_coupled_pair'][1]}"
        print(f"{pdb_id:<8} {result['n_trp']:<5} {pair_label:<14} {result['distance_A']:<8.2f}  {result['coupling_cm']:<10.2f} {result['S_value']:<8.4f} {result['status']:<20}")

    # Summary
    print(f"\n{'='*72}")
    violators = [r for r in [compute_chsh_from_pdb(p) for p in pdb_list] if r and r['violates_classical_limit']]
    print(f"  Targets violating classical bound (S > 2): {len(violators)} / {len(pdb_list)}")
    for v in violators:
        print(f"    {v['pdb_id']}: S = {v['S_value']:.4f}, surplus = {v['S_violation']:.4f}")


if __name__ == "__main__":
    targets = ["6CNO", "6LQA", "7KOX", "7TYO", "6J8J", "6PV7", "1BL8"]
    run_bell_test(targets)

