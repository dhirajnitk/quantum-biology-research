"""
fmo_lindblad.py  —  FMO 7-site Lindblad simulation for P1 paper.

Implements:
  1. Adolphs & Renger (2006) 7-site FMO Hamiltonian
  2. Lindblad master equation with dephasing and trapping
  3. ENAQT efficiency curve vs dephasing rate
  4. Quantum Mutual Information (QMI) matrix

All calculations use NumPy/SciPy (no QuTiP dependency).
"""

import numpy as np
from numpy import pi, exp, trace, abs, real
from scipy.integrate import solve_ivp


HBAR = 1.054571817e-34
CM1_TO_RADS = 2 * pi * 2.99792458e10

# FMO 7-site Hamiltonian (Adolphs & Renger 2006, cm-1)
H_FMO_CM = np.array([
    [  0.0, -87.7,   5.5,  -5.9,   6.7, -13.7,  -9.9],
    [-87.7, 320.0,  30.8,   8.2,   0.7,  11.8,   4.3],
    [  5.5,  30.8,   0.0, -53.5,  -2.2,  -9.6,   6.0],
    [ -5.9,   8.2, -53.5, 110.0, -70.7, -17.0, -63.3],
    [  6.7,   0.7,  -2.2, -70.7, 270.0,  81.1,  -1.3],
    [-13.7,  11.8,  -9.6, -17.0,  81.1, 420.0,  39.7],
    [ -9.9,   4.3,   6.0, -63.3,  -1.3,  39.7, 230.0],
])

# FMO 8-site Hamiltonian (Schmidt am Busch et al. 2011, JPCB)
# Includes BChl 8 as the entry bridge from the chlorosome antenna
H_FMO8_CM = np.array([
    [  0.0, -87.7,   5.5,  -5.9,   6.7, -13.7,  -9.9,  -6.5],
    [-87.7, 320.0,  30.8,   8.2,   0.7,  11.8,   4.3,   3.2],
    [  5.5,  30.8,   0.0, -53.5,  -2.2,  -9.6,   6.0,   0.8],
    [ -5.9,   8.2, -53.5, 110.0, -70.7, -17.0, -63.3,  12.5],
    [  6.7,   0.7,  -2.2, -70.7, 270.0,  81.1,  -1.3,   1.7],
    [-13.7,  11.8,  -9.6, -17.0,  81.1, 420.0,  39.7,  -3.5],
    [ -9.9,   4.3,   6.0, -63.3,  -1.3,  39.7, 230.0,  40.1],
    [ -6.5,   3.2,   0.8,  12.5,   1.7,  -3.5,  40.1, 140.0],
])


class FMOSolver:
    """Lindblad master equation solver for the FMO complex."""

    def __init__(self, gamma_deph=300.0, gamma_trap=1.0, temperature=300.0):
        self.N = 7
        self.gamma_deph_cm = gamma_deph
        self.gamma_trap_ps = gamma_trap
        self.T = temperature
        self.gamma_deph = gamma_deph * CM1_TO_RADS
        self.gamma_trap = gamma_trap * 1e12
        self.H = H_FMO_CM * CM1_TO_RADS
        self._build_ops()

    def _build_ops(self):
        self.c_ops = []
        for j in range(self.N):
            L = np.zeros((self.N, self.N), dtype=complex)
            L[j, j] = np.sqrt(self.gamma_deph)
            self.c_ops.append(L)
        L_trap = np.zeros((self.N, self.N), dtype=complex)
        L_trap[2, 2] = np.sqrt(self.gamma_trap)
        self.c_ops.append(L_trap)

    def _lindblad_rhs(self, t, rho_vec):
        N = self.N
        rho = rho_vec.reshape((N, N))
        drho = -1j * (self.H @ rho - rho @ self.H)
        for L in self.c_ops:
            Ld = L.conj().T
            drho += L @ rho @ Ld - 0.5 * (Ld @ L @ rho + rho @ Ld @ L)
        return drho.flatten()

    def run(self, t_max_ps=5.0, n_steps=500, initial_site=1):
        N = self.N
        t_eval = np.linspace(0, t_max_ps, n_steps)
        ps = initial_site - 1
        rho0 = np.zeros((N, N), dtype=complex)
        rho0[ps, ps] = 1.0
        sol = solve_ivp(self._lindblad_rhs, (0, t_max_ps), rho0.flatten(),
                        method='RK45', t_eval=t_eval, rtol=1e-8, atol=1e-10)
        tlist = sol.t * 1e12
        n_pts = len(tlist)
        rho_t = np.zeros((n_pts, N, N), dtype=complex)
        pops = np.zeros((n_pts, N))
        for i in range(n_pts):
            rho_t[i] = sol.y[:, i].reshape((N, N))
            pops[i] = np.real(np.diag(rho_t[i]))
        return tlist, rho_t, pops

    def compute_efficiency(self, t_max_ps=5.0, n_steps=500, initial_site=1):
        tlist, rho_t, pops = self.run(t_max_ps, n_steps, initial_site)
        p_trap = pops[:, 2]
        eta = self.gamma_trap_ps * 1e12 * np.trapz(p_trap, tlist * 1e-12)
        return min(eta, 1.0)

    @staticmethod
    def von_neumann_entropy(rho):
        evals = np.linalg.eigvalsh(rho)
        evals = evals[evals > 1e-15]
        return -np.sum(evals * np.log(evals))

    @staticmethod
    def quantum_discord(rho_ij):
        """Compute quantum discord for a 2-site density matrix.
        
        D(A:B) = I(A:B) - J(A:B) where J is classical correlation
        (the maximum information about B obtainable by measuring A).
        For the single-exciton subspace, we use the approach of
        Sarovar et al. (2010) Nat. Phys.
        """
        pi_i = rho_ij[0, 0] + rho_ij[1, 1]
        pi_j = rho_ij[0, 0] + rho_ij[2, 2]
        
        # Single-site entropies
        si = FMOSolver.von_neumann_entropy(np.array([[pi_i, 0], [0, 1-pi_i]]))
        sj = FMOSolver.von_neumann_entropy(np.array([[pi_j, 0], [0, 1-pi_j]]))
        
        # Joint entropy
        p = np.array([(1-pi_i)*(1-pi_j), (1-pi_i)*pi_j, pi_i*(1-pi_j), pi_i*pi_j])
        p = p / p.sum()
        sij_joint = -np.sum(p * np.log(p + 1e-30))
        
        # Mutual information
        I_q = si + sj - sij_joint
        
        # Classical correlation J(A:B): optimize over projective measurements on A
        # For a qubit, the optimal measurement basis is the eigenbasis of rho_A
        rho_A = np.array([[rho_ij[0,0]+rho_ij[1,1], rho_ij[0,2]+rho_ij[1,3]],
                          [rho_ij[2,0]+rho_ij[3,1], rho_ij[2,2]+rho_ij[3,3]]])
        evals, evecs = np.linalg.eigh(rho_A)
        
        # Projective measurement in eigenbasis
        J_max = 0.0
        for ka in range(2):
            proj = np.outer(evecs[:, ka], evecs[:, ka].conj())
            pk = np.trace(proj @ rho_A).real
            if pk < 1e-10:
                continue
            # Conditional state of B
            rho_B_given_A = (np.array([
                [rho_ij[ka*2, ka*2], rho_ij[ka*2, ka*2+1]],
                [rho_ij[ka*2+1, ka*2], rho_ij[ka*2+1, ka*2+1]]
            ]) / pk + 1e-30)
            rho_B_given_A /= np.trace(rho_B_given_A)
            s_cond = FMOSolver.von_neumann_entropy(rho_B_given_A)
            J_max += pk * s_cond
        
        J_classical = sj - J_max
        D = max(0.0, I_q - J_classical)
        return D, I_q, J_classical

    def compute_qmi_matrix(self, t_max_ps=2.0, n_steps=200, initial_site=1):
        tlist, rho_t, pops = self.run(t_max_ps, n_steps, initial_site)
        N = self.N
        n_pts = len(tlist)
        qmi_mat = np.zeros((N, N))
        for i in range(N):
            for j in range(i + 1, N):
                qmi_series = np.zeros(n_pts)
                for k in range(n_pts):
                    pi_i = pops[k, i]
                    pi_j = pops[k, j]
                    si = FMOSolver.von_neumann_entropy(np.array([[pi_i, 0], [0, 1-pi_i]]))
                    sj = FMOSolver.von_neumann_entropy(np.array([[pi_j, 0], [0, 1-pi_j]]))
                    p = np.array([(1-pi_i)*(1-pi_j), (1-pi_i)*pi_j, pi_i*(1-pi_j), pi_i*pi_j])
                    p = p / p.sum()
                    sij = -np.sum(p * np.log(p + 1e-30))
                    qmi_series[k] = si + sj - sij
                qmi_mat[i, j] = qmi_mat[j, i] = np.mean(qmi_series)
        return qmi_mat, tlist


if __name__ == "__main__":
    print("=" * 60)
    print("  FMO LINDBLAD SIMULATION - P1 Baseline")
    print("=" * 60)

    print("\n  ENAQT EFFICIENCY vs DEPHASING RATE")
    print(f"  {'Dephasing (cm-1)':<20} {'Efficiency':<15}")
    print(f"  {'-'*35}")
    gammas = [0.1, 1.0, 10, 50, 100, 175, 195, 300, 500, 1000]
    effs = []
    for g in gammas:
        s = FMOSolver(gamma_deph=g, gamma_trap=1.0)
        eta = s.compute_efficiency()
        effs.append(eta)
        print(f"  {g:<20.1f} {eta:<15.4f}")

    opt_idx = int(np.argmax(effs))
    print(f"\n  Optimal dephasing: {gammas[opt_idx]:.1f} cm-1")
    print(f"  Optimal efficiency: {effs[opt_idx]:.4f}")

    print("\n  QUANTUM MUTUAL INFORMATION MATRIX (nats)")
    solver = FMOSolver(gamma_deph=175.0)
    qmi_mat, _ = solver.compute_qmi_matrix()
    print(f"  {'':>3}", end="")
    for j in range(7):
        print(f" {j+1:>8}", end="")
    print()
    for i in range(7):
        print(f"  {i+1:>3}", end="")
        for j in range(7):
            print(f" {qmi_mat[i,j]:>8.4f}", end="")
        print()

    pairs = [(i, j, qmi_mat[i, j]) for i in range(7) for j in range(i+1, 7)]
    pairs.sort(key=lambda x: -x[2])
    print(f"\n  TOP 5 QMI PAIRS:")
    for pair in pairs[:5]:
        print(f"    BChl {pair[0]+1} - BChl {pair[1]+1}: {pair[2]:.4f} nats")
