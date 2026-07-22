"""
lindblad_solver.py  —  FMO 7-site Lindblad solver for ENAQT and QMI.

Implements:
  - Adolphs & Renger (2006) 7-site FMO Hamiltonian
  - Liouville-space Lindblad master equation via matrix exponentiation
  - ENAQT efficiency curve (correct: low at gamma=0 due to Anderson localization)
  - Quantum Mutual Information matrix

Trapping is handled via a non-Hermitian Hamiltonian term.
Pure dephasing via standard Lindblad |k><k| operators.

References:
  Adolphs & Renger (2006) Biophys J 91, 2778-2797
  Mohseni et al. (2008) JCP 129, 174106
  Ishizaki & Fleming (2009) PNAS 106, 17255
"""

import numpy as np
from scipy.linalg import expm


class FmoLindbladSolver:
    """FMO 7-site Lindblad solver with Liouville-space propagation.

    Uses non-Hermitian Hamiltonian for trapping and Lindblad dephasing.
    """

    def __init__(self, disorder_std=0.0):
        # Adolphs & Renger (2006) 7-site Hamiltonian, cm^-1
        self.H_cm = np.array([
            [215, -104.1,   5.5,  -5.9,   4.3, -15.1,  -6.7],
            [-104.1, 220,   28.2,   7.1,  -5.3,   5.6,   4.4],
            [  5.5,  28.2,   0.0, -60.8,  -1.5,  -0.2,   0.9],
            [ -5.9,   7.1, -60.8, 125.0, -73.4,  -2.0, -61.2],
            [  4.3,  -5.3,  -1.5, -73.4, 450.0,  43.0, -13.8],
            [-15.1,   5.6,  -0.2,  -2.0,  43.0, 310.0,  34.6],
            [ -6.7,   4.4,   0.9, -61.2, -13.8,  34.6, 230.0],
        ], dtype=float)

        # Add static disorder to site energies (Anderson localization)
        # Disorder suppresses transport at gamma=0, creating ENAQT peak
        if disorder_std > 0:
            rng = np.random.RandomState(42)
            disorder = rng.normal(0, disorder_std, 7)
            for i in range(7):
                self.H_cm[i, i] += disorder[i]

        # Alias for API compatibility
        self.H_base = self.H_cm

        # Convert to ps^-1 (hbar=1): 1 cm^-1 = 0.0188 ps^-1
        self.H = self.H_cm * 0.0188
        self.n_sites = 7
        self.disorder_std = disorder_std
        # Trapping rate at site 3 (index 2)
        self.k_trap = 1.0  # ps^-1

    def _build_lindbladian(self, dephasing_rate):
        """Construct Liouville-space Lindblad superoperator.

        L = -i[H_eff, .] + sum_k gamma_k (L_k . L_k^dag - 1/2 {L_k^dag L_k, .})

        H_eff = H - i * k_trap/2 * |2><2|  (non-Hermitian trapping)
        L_k = |k><k|  (pure dephasing at each site)
        """
        n = self.n_sites
        dim = n ** 2
        L_super = np.zeros((dim, dim), dtype=complex)
        H = self.H.copy()

        # Non-Hermitian trapping: H_eff = H - i * k_trap/2 * |2><2|
        H = H.astype(complex)
        H[2, 2] -= 1j * self.k_trap / 2.0

        # Unitary evolution with non-Hermitian H: -i(H rho - rho H†)
        # In Liouville space: vec(-iH rho) = -i(I ⊗ H) vec(rho)
        #                    vec(i rho H†) = i(H†ᵀ ⊗ I) vec(rho) = i(H̅ ⊗ I) vec(rho)
        # where H̅ is the complex conjugate of H (NOT the transpose)
        L_super = -1j * (np.kron(np.eye(n), H) - np.kron(H.conj(), np.eye(n)))

        # Pure dephasing at each site: L_k = |k><k|
        for k in range(n):
            P_k = np.zeros((n, n))
            P_k[k, k] = 1.0
            # L_k rho L_k^dag -> (L_k^T ⊗ L_k) vec(rho) = (P_k ⊗ P_k) vec(rho)
            # -1/2 {L_k^dag L_k, rho} -> -1/2 (I ⊗ P_k + P_k^T ⊗ I) vec(rho)
            deph_super = np.kron(P_k, P_k) - 0.5 * (np.kron(np.eye(n), P_k) + np.kron(P_k.T, np.eye(n)))
            L_super += dephasing_rate * deph_super

        return L_super

    def run_time_evolution(self, dephasing_rate, t_max=10.0, dt=0.01):
        """Run Liouville-space Lindblad evolution.

        Efficiency computed from population lost through the trap:
        eta = k_trap * integral of rho[2,2] over time

        Returns:
            efficiency: fraction of excitons trapped at reaction center
            density_history: list of density matrices at each time step
        """
        L = self._build_lindbladian(dephasing_rate)
        n = self.n_sites
        t_steps = int(t_max / dt)

        rho_0 = np.zeros((n, n), dtype=complex)
        rho_0[0, 0] = 1.0
        vec = rho_0.flatten().reshape(-1, 1)

        propagator = expm(L * dt)
        trapped = 0.0
        history = []

        for _ in range(t_steps):
            vec = propagator @ vec
            rho = vec.reshape((n, n))
            history.append(rho.copy())
            trapped += self.k_trap * np.real(rho[2, 2]) * dt

        return min(trapped, 1.0), history

    def compute_qmi_matrix(self, rho):
        """Quantum Mutual Information matrix for all site pairs."""
        n = self.n_sites
        pops = np.real(np.diag(rho))
        pops = np.clip(pops, 1e-12, 1.0)
        s_single = -pops * np.log2(pops)

        qmi = np.zeros((n, n))
        for i in range(n):
            qmi[i, i] = s_single[i]
            for j in range(i + 1, n):
                rho_sub = np.array([[pops[i], rho[i, j]],
                                    [rho[j, i], pops[j]]])
                ev = np.clip(np.linalg.eigvalsh(rho_sub), 1e-12, 1.0)
                s_joint = -np.sum(ev * np.log2(ev))
                val = max(0.0, s_single[i] + s_single[j] - s_joint)
                qmi[i, j] = qmi[j, i] = val
        return qmi

    # API compatibility alias
    calculate_quantum_mutual_information = compute_qmi_matrix

    def enaqt_curve(self, gamma_range=None):
        """Compute ENAQT efficiency vs dephasing rate (cm^-1)."""
        if gamma_range is None:
            gamma_range = [0, 10, 50, 100, 175, 195, 300, 500, 1000]
        effs = []
        for g_cm in gamma_range:
            g_ps = g_cm * 0.0188  # cm^-1 -> ps^-1
            eff, _ = self.run_time_evolution(dephasing_rate=g_ps)
            effs.append(eff)
        return gamma_range, effs


if __name__ == "__main__":
    solver = FmoLindbladSolver()
    print("  FMO ENAQT EFFICIENCY vs DEPHASING RATE")
    print(f"  {'Gamma (cm-1)':<16} {'Gamma (ps-1)':<16} {'Efficiency':<15}")
    print(f"  {'-'*47}")
    gammas, effs = solver.enaqt_curve()
    for g_cm, eff in zip(gammas, effs):
        g_ps = g_cm * 0.0188
        print(f"  {g_cm:<16.0f} {g_ps:<16.3f} {eff:<15.4f}")

    opt_idx = int(np.argmax(effs))
    print(f"\n  Optimal dephasing: {gammas[opt_idx]:.0f} cm-1 ({gammas[opt_idx]*0.0188:.2f} ps-1)")
    print(f"  Peak efficiency: {effs[opt_idx]*100:.1f}%")

    # QMI at ENAQT optimal
    eff, hist = solver.run_time_evolution(dephasing_rate=175 * 0.0188)
    qmi = solver.compute_qmi_matrix(hist[-1])
    print(f"\n  QMI MATRIX (bits) at gamma=175 cm-1, t=10ps")
    for i in range(7):
        row = '  '.join(f"{qmi[i,j]:.4f}" for j in range(7))
        print(f"    {row}")
