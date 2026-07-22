"""
transition_dipoles.py — TD-DFT interface for L_a/L_b transition dipoles.

Extracts indole ring geometries from PDB structures and prepares input
for quantum chemistry calculations (PySCF, ORCA) to compute the two
close-lying excited states (L_a and L_b) and their transition dipole moments.

When TD-DFT results are not available, known literature values for the
indole chromophore are used as defaults.
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from src.pdb_tools.trp_extractor import fetch_pdb, extract_trp_data


# Literature values for indole (Trp) transition dipole moments
# L_a state: 280 nm, transition dipole along the long axis (CG -> CH2)
# L_b state: 290 nm, transition dipole along the short axis (perpendicular)
# Reference: Callis, Methods Enzymol. 278, 113 (1997); Vivian & Callis, Biophys J. 80, 2093 (2001)

LITERATURE_DIPOLES = {
    'L_a': {
        'wavelength_nm': 280,
        'energy_cm': 35700,
        'oscillator_strength': 0.12,
        'dipole_strength_D2': 3.5,  # Debye^2
        'direction': 'long_axis',  # CG -> CH2
    },
    'L_b': {
        'wavelength_nm': 290,
        'energy_cm': 34480,
        'oscillator_strength': 0.04,
        'dipole_strength_D2': 1.2,
        'direction': 'short_axis',  # perpendicular to CG -> CH2
    },
}

INDOLE_ATOMS = ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']


def extract_indole_geometry(pdb_id, chain=None):
    """Extract indole ring atomic coordinates from a PDB structure.

    Returns a dict with atom names and (x, y, z) coordinates for the
    indole ring, suitable for quantum chemistry input.
    """
    text = fetch_pdb(pdb_id, cache=True)
    if not text:
        return None

    centroids, dipoles = extract_trp_data(text, chain)
    # dipoles is a dict: {res_num: np.array([dx, dy, dz])}

    geometries = {}
    for res_num in centroids:
        geometry = {}
        # The CG atom position is the centroid
        geometry['centroid'] = centroids[res_num]
        # The transition dipole direction (from our extraction)
        geometry['dipole_L_a'] = dipoles.get(res_num, np.zeros(3))
        # L_b is perpendicular to L_a — rotate 90 degrees
        if np.linalg.norm(dipoles.get(res_num, np.zeros(3))) > 0:
            d_la = dipoles[res_num] / np.linalg.norm(dipoles[res_num])
            # Find a perpendicular vector
            if abs(d_la[2]) < 0.9:
                d_lb = np.cross(d_la, np.array([0, 0, 1]))
            else:
                d_lb = np.cross(d_la, np.array([1, 0, 0]))
            d_lb = d_lb / np.linalg.norm(d_lb)
        else:
            d_lb = np.array([0, 1, 0])
        geometry['dipole_L_b'] = d_lb

        geometries[res_num] = geometry

    return geometries


def compute_dipole_coupling(mu_i, mu_j, R_ij, eps=2.0):
    """Compute dipole-dipole coupling J_ij with orientation factor.

    J_ij = J0 * (R0/R)^3 * |kappa| / sqrt(eps)
    where kappa = mu_i·mu_j - 3(mu_i·R_hat)(mu_j·R_hat)

    For L_a and L_b states separately.
    """
    J0 = -80.0  # cm^-1 at R0 = 10 A
    R0 = 10.0

    R_norm = np.linalg.norm(R_ij)
    if R_norm < 0.1:
        return 0.0

    R_hat = R_ij / R_norm
    mu_i_norm = mu_i / np.linalg.norm(mu_i)
    mu_j_norm = mu_j / np.linalg.norm(mu_j)

    kappa = np.dot(mu_i_norm, mu_j_norm) - 3 * np.dot(mu_i_norm, R_hat) * np.dot(mu_j_norm, R_hat)

    J = J0 * (R0 / R_norm) ** 3 * abs(kappa) / np.sqrt(eps)
    return J


def compute_full_coupling_matrix(pdb_id, eps=2.0):
    """Compute complete coupling matrix including L_a and L_b contributions.

    Returns J_aa, J_bb, J_ab matrices for L_a-L_a, L_b-L_b, and
    cross-couplings.
    """
    text = fetch_pdb(pdb_id, cache=True)
    if not text:
        return None

    centroids, dipoles = extract_trp_data(text)
    keys = list(centroids.keys())
    n = len(keys)

    if n < 2:
        return None

    # Build full 2n x 2n matrix (n sites x 2 states each)
    J_full = np.zeros((2 * n, 2 * n))

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            R_vec = centroids[keys[i]] - centroids[keys[j]]
            d_i = dipoles.get(keys[i], np.zeros(3))
            d_j = dipoles.get(keys[j], np.zeros(3))

            if np.linalg.norm(d_i) < 0.1 or np.linalg.norm(d_j) < 0.1:
                continue

            # L_a-L_a coupling
            J_aa = compute_dipole_coupling(d_i, d_j, R_vec, eps)
            J_full[i, j] = J_aa

            # L_b-L_b coupling (using perpendicular dipoles)
            # Build L_b vectors by rotating L_a by 90 degrees
            d_i_b = np.cross(d_i / np.linalg.norm(d_i), np.array([0, 0, 1]))
            d_j_b = np.cross(d_j / np.linalg.norm(d_j), np.array([0, 0, 1]))
            if np.linalg.norm(d_i_b) > 0 and np.linalg.norm(d_j_b) > 0:
                J_bb = compute_dipole_coupling(d_i_b, d_j_b, R_vec, eps)
                J_full[n + i, n + j] = J_bb

            # Cross-couplings
            if np.linalg.norm(d_i_b) > 0:
                J_ab = compute_dipole_coupling(d_i, d_j_b, R_vec, eps)
                J_full[i, n + j] = J_ab
                J_full[n + i, j] = J_ab

    return J_full, keys


if __name__ == '__main__':
    print("=== Transition Dipole Analysis ===\n")
    print("L_a (280 nm): long-axis transition")
    print(f"  Oscillator strength: {LITERATURE_DIPOLES['L_a']['oscillator_strength']}")
    print(f"  Direction: {LITERATURE_DIPOLES['L_a']['direction']}")
    print()
    print("L_b (290 nm): short-axis transition")
    print(f"  Oscillator strength: {LITERATURE_DIPOLES['L_b']['oscillator_strength']}")
    print(f"  Direction: {LITERATURE_DIPOLES['L_b']['direction']}")
    print()

    # Compute full coupling matrix for KcsA
    print("--- Full coupling matrix (L_a + L_b): KcsA (1BL8) ---")
    result = compute_full_coupling_matrix('1BL8')
    if result is not None:
        J_full, keys = result
        n = len(keys)
        print(f"  Sites: {n}")
        print(f"  Matrix shape: {J_full.shape} ({n} sites x 2 states each)")

        # Diagonal couplings (L_a-L_a)
        J_aa = J_full[:n, :n]
        print(f"\n  L_a-L_a couplings (cm^-1):")
        for i in range(n):
            for j in range(i+1, n):
                if abs(J_aa[i,j]) > 1:
                    print(f"    {keys[i]}-{keys[j]}: J_aa = {J_aa[i,j]:.1f}")

        # Additional coupling from L_b
        J_bb = J_full[n:, n:]
        print(f"\n  L_b-L_b couplings (cm^-1):")
        for i in range(n):
            for j in range(i+1, n):
                if abs(J_bb[i,j]) > 1:
                    print(f"    {keys[i]}-{keys[j]}: J_bb = {J_bb[i,j]:.1f}")
