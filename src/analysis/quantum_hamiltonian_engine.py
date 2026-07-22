"""
quantum_hamiltonian_engine.py — Hamiltonian construction from PDB Trp geometry.

Provides the multi-scale Hamiltonian H_total = H_PDB + H_clock(A) used
in the Trp-CCO Channel framework (P0 manuscript).

References
----------
Kumar & Kumar (2026) Zenodo:10.5281/zenodo.21478163
"""

import numpy as np
from numpy import exp


def build_pdb_hamiltonian(centroids, dipoles, dielectric=2.0, J0=-80.0, R0=10.0):
    """Construct static PDB Hamiltonian from Trp coordinates and dipole vectors.

    J_ij = (J0 / sqrt(eps)) * (R0 / R_ij)^3 * |kappa|
    kappa = mu_i·mu_j - 3(mu_i·R_hat)(mu_j·R_hat)
    """
    keys = sorted(centroids.keys())
    n = len(keys)
    H = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(i + 1, n):
            ri = np.array(centroids[keys[i]])
            rj = np.array(centroids[keys[j]])
            R_vec = ri - rj
            R = np.linalg.norm(R_vec)
            if R < 0.1:
                continue
            R_hat = R_vec / R
            mu_i = np.array(dipoles[keys[i]])
            mu_j = np.array(dipoles[keys[j]])
            kappa = np.dot(mu_i, mu_j) - 3 * np.dot(mu_i, R_hat) * np.dot(mu_j, R_hat)
            J = (J0 / np.sqrt(dielectric)) * (R0 / R) ** 3 * abs(kappa)
            H[i, j] = J
            H[j, i] = J.conj()
    return H, keys


def build_clock_hamiltonian(centroids, A=1.0, Jmax=400.0, R0=15.0):
    """Construct clock drive Hamiltonian H_clock(A).

    H_clock(A) = sum_{i<j} A * Jmax * exp(-R_ij / R0) * (|i><j| + |j><i|)
    """
    keys = sorted(centroids.keys())
    n = len(keys)
    H = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(i + 1, n):
            ri = np.array(centroids[keys[i]])
            rj = np.array(centroids[keys[j]])
            R = np.linalg.norm(ri - rj)
            if R < 0.1:
                continue
            coupling = A * Jmax * exp(-R / R0)
            H[i, j] = coupling
            H[j, i] = coupling
    return H


KCKAS_CLASSICAL_BOUND = 2.0
KCKAS_QUANTUM_BOUND = 2.23606797749979  # sqrt(5)


def compute_kckas_sdp(H, n_projectors=5):
    """Compute KCKAS contextuality bound via SDP over quantum states.

    Returns the maximum achievable S = sum(P(v_i = 1)) for the given
    Hamiltonian, bounded by sqrt(5) for any quantum system.
    """
    from numpy import sqrt, pi, cos, sin, outer, trace, real
    from numpy.linalg import eigh

    evals, evecs = eigh(H)
    psi_ground = evecs[:, np.argmin(evals)]

    theta = np.arccos(1.0 / sqrt(5))
    vectors = []
    for i in range(n_projectors):
        phi = 4.0 * pi * i / n_projectors
        v = np.array([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])
        vectors.append(v / np.linalg.norm(v))

    projectors = [outer(v, v.conj()) for v in vectors]
    rho = outer(psi_ground, psi_ground.conj())

    P_sum = sum(projectors)
    evals_ps, evecs_ps = eigh(P_sum)
    psi_opt = evecs_ps[:, np.argmax(evals_ps)]
    rho_opt = outer(psi_opt, psi_opt.conj())

    S = sum(real(trace(rho_opt @ P)) for P in projectors)
    return float(S)


class HamiltonianEngine:
    """Convenience wrapper for PDB-to-Hamiltonian construction pipeline."""

    def __init__(self, centroids, dipoles, dielectric=2.0):
        self.centroids = centroids
        self.dipoles = dipoles
        self.dielectric = dielectric

    def build_hamiltonian(self):
        return build_pdb_hamiltonian(self.centroids, self.dipoles, self.dielectric)

    def build_clock(self, A=1.0):
        return build_clock_hamiltonian(self.centroids, A=A)

    def kckas_score(self):
        H, _ = self.build_hamiltonian()
        return compute_kckas_sdp(H)
