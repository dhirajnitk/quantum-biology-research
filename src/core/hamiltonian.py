"""
hamiltonian.py  —  Formal multi-scale Hamiltonian for PDB-grounded Trp networks.

H_total = H_PDB + H_clock(A)

  H_PDB:   Static structural baseline from PDB coordinates.
           J_ij = (J0 / sqrt(eps)) * (R0 / R_ij)^3 * kappa  (Dexter dipole-dipole)
           where kappa = mu_i·mu_j - 3(mu_i·R_hat)(mu_j·R_hat) is the
           dipole-dipole orientation factor (ranges from -2 to 2).

  H_clock: Macro-scale oscillatory drive from Layer 3 (gamma rhythms).
           J_clock(R) = A * J_max * exp(-R / R_c)  (exponential phase pump)

References
----------
- Adolphs & Renger (2006) Biophys J
- Firmenich et al. (2026) bioRxiv
- Valeur & Berberan-Santos (2012) Molecular Fluorescence (dipole orientation)
"""

import numpy as np
from numpy import exp, sqrt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))


# ── Default physical parameters ─────────────────────────────────────
J0 = -80.0          # cm-1, reference coupling at R0
R0 = 10.0           # A, reference distance
J_MAX = 400.0       # cm-1, maximum clock-driven coupling
R_CLOCK = 15.0      # A, clock coupling decay length
DISORDER_STD = 30.0 # cm-1, static site energy disorder
SITE_SPACING = 80.0 # cm-1, mean site energy separation


def calculate_orientation_factor(mu_i, mu_j, R_ij):
    """Dipole-dipole orientation factor kappa between two Trp rings.

    kappa = mu_i_hat . mu_j_hat - 3 (mu_i_hat . R_hat)(mu_j_hat . R_hat)

    |kappa| ranges from 0 to 2.
    kappa = 0 when dipoles are orthogonal or at magic angle.
    |kappa| = 2 for collinear head-to-tail alignment.

    Parameters
    ----------
    mu_i, mu_j : ndarray
        Transition dipole unit vectors for rings i and j.
    R_ij : ndarray
        Displacement vector between ring centres (Angstroms).

    Returns
    -------
    kappa : float
        Orientation factor in [-2, 2].
    """
    r_hat = R_ij / max(np.linalg.norm(R_ij), 1e-30)
    mu_i_u = mu_i / max(np.linalg.norm(mu_i), 1e-30)
    mu_j_u = mu_j / max(np.linalg.norm(mu_j), 1e-30)
    kappa = np.dot(mu_i_u, mu_j_u) - 3 * np.dot(mu_i_u, r_hat) * np.dot(mu_j_u, r_hat)
    return kappa


def build_pdb_hamiltonian(xyz_coords, dielectric=2.0, site_energies=None,
                          dipole_vectors=None):
    """Build H_PDB from PDB-derived Trp coordinates.

    Parameters
    ----------
    xyz_coords : ndarray
        (N, 3) array of Trp ring centroid coordinates in Angstroms.
    dielectric : float
        Local dielectric constant (default 2.0 for hydrophobic core).
    site_energies : ndarray or None
        Site energies in cm-1. If None, uses SITE_SPACING ramp + disorder.

    Returns
    -------
    H_cm : ndarray
        (N, N) Hamiltonian matrix in cm-1.
    distances : ndarray
        (N, N) pair-distance matrix in Angstroms.
    """
    n = len(xyz_coords)
    H = np.zeros((n, n), dtype=complex)
    D = np.zeros((n, n))

    # Build distance matrix
    for i in range(n):
        for j in range(n):
            D[i, j] = np.linalg.norm(xyz_coords[i] - xyz_coords[j])

    # Site energies
    if site_energies is None:
        rng = np.random.RandomState(42)
        site_energies = np.arange(n) * SITE_SPACING + rng.normal(0, DISORDER_STD, n)

    for i in range(n):
        H[i, i] = site_energies[i]

    # Off-diagonal couplings (Dexter-type with dipole orientation)
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] > 0.1:
                J_ij = (J0 / sqrt(dielectric)) * (R0 / D[i, j]) ** 3

                # Apply dipole orientation factor if available
                if dipole_vectors is not None and i < len(dipole_vectors) and j < len(dipole_vectors):
                    R_vec = xyz_coords[j] - xyz_coords[i]
                    kappa = calculate_orientation_factor(
                        dipole_vectors[i], dipole_vectors[j], R_vec)
                    J_ij *= abs(kappa)

                H[i, j] = J_ij
                H[j, i] = J_ij

    return H, D


def build_clock_hamiltonian(distances, amplitude=1.0):
    """Build H_clock(A) — macro-scale clock drive.

    Adds coherent off-diagonal coupling that overcomes site energy
    disorder and aligns fluctuating Trp dipole phases.

    Parameters
    ----------
    distances : ndarray
        (N, N) pair-distance matrix in Angstroms.
    amplitude : float
        Clock amplitude A in [0, 1].

    Returns
    -------
    H_clock : ndarray
        (N, N) clock-driven coupling matrix in cm-1.
    """
    n = len(distances)
    H = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(i + 1, n):
            if distances[i, j] > 0.1:
                J_clock = amplitude * J_MAX * exp(-distances[i, j] / R_CLOCK)
                H[i, j] = J_clock
                H[j, i] = J_clock
    return H


def build_multiscale_hamiltonian(xyz_coords, clock_amplitude=0.0,
                                  dielectric=2.0, site_energies=None,
                                  dipole_vectors=None):
    """Build the complete multi-scale Hamiltonian.

    H_total = H_PDB + H_clock(A)

    Parameters
    ----------
    xyz_coords : ndarray
        (N, 3) array of Trp centroid coordinates in Angstroms.
    clock_amplitude : float
        A in [0, 1]. 0 = static crystal, 1 = fully synchronised.
    dielectric : float
        Local dielectric constant.
    site_energies : ndarray or None

    Returns
    -------
    H_total_cm : ndarray
        (N, N) total Hamiltonian in cm-1.
    distances : ndarray
        (N, N) distance matrix in Angstroms.
    """
    H_pdb, D = build_pdb_hamiltonian(xyz_coords, dielectric, site_energies, dipole_vectors)
    H_clock = build_clock_hamiltonian(D, clock_amplitude)
    return H_pdb + H_clock, D


def coherence_factor(H_cm):
    """Compute coherence factor F_coh from Hamiltonian matrix.

    F_coh = sum(|H_ij|) / (sum(|H_ij|) + sum(|H_ii|))

    Ranges from 0 (fully classical) to ~1 (fully coherent).
    """
    n = H_cm.shape[0]
    off_sum = 0.0
    diag_sum = 0.0
    for i in range(n):
        diag_sum += abs(H_cm[i, i])
        for j in range(n):
            if i != j:
                off_sum += abs(H_cm[i, j])
    if diag_sum + off_sum < 1e-30:
        return 0.0
    return off_sum / (off_sum + diag_sum)


def kckas_s(H_cm):
    """KCKAS coherence score S from Hamiltonian.

    S = F_coh * sqrt(5)
    NOTE: This is a geometric proxy based on the coherence factor.
    The true KCKAS contextuality sum requires optimal state preparation
    (see analysis/pdb_contextuality.py for the full implementation).

    Maximum quantum value: sqrt(5) = 2.2361
    Classical bound: S <= 2.0
    """
    F = coherence_factor(H_cm)
    return F * np.sqrt(5)
