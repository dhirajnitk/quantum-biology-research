"""
quantum_darwinism.py  —  Quantum Darwinism analysis for FMO complex.

Implements Zurek's Quantum Darwinism framework for the FMO 7-site complex:
  1. Partial Information Curves (PIC) — I(rho_sys : rho_frag) vs fragment size
  2. Redundancy R(delta) — number of independent fragments with I >= delta
  3. Pointer states — states that survive decoherence (einselection)
  4. Spectrum Broadcast Structure (SBS) — classical information capacity

Theory:
  Zurek (2003) RMP 75, 715 — Decoherence, einselection, QD
  Horodecki et al. (2015) Nature Comms — SBS framework
"""

import numpy as np
from numpy import log2, clip, real, diag, linalg, zeros, ix_
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.core.lindblad_solver import FmoLindbladSolver


class QuantumDarwinism:
    """Complete Quantum Darwinism analysis for the FMO complex."""

    def __init__(self, density_history, system_site=0):
        self.history = density_history
        self.n_sites = density_history[0].shape[0] if density_history else 7
        self.sys_idx = system_site  # system = input site (default site 1)
        self.env_indices = [i for i in range(self.n_sites) if i != self.sys_idx]

    def von_neumann_entropy(self, rho):
        ev = linalg.eigvalsh(rho)
        ev = clip(ev, 1e-15, 1.0)
        return -np.sum(ev * log2(ev))

    def mutual_information(self, rho_sys, rho_env, rho_joint):
        """I(sys:env) = S(rho_sys) + S(rho_env) - S(rho_joint)."""
        return max(0.0, self.von_neumann_entropy(rho_sys) +
                   self.von_neumann_entropy(rho_env) -
                   self.von_neumann_entropy(rho_joint))

    def compute_redundancy_curve(self, time_index):
        """Partial Information Curve: I(system : fragment) vs fragment size.

        System = site `sys_idx` (default site 1, the input gateway).
        Fragment = subsets of remaining sites.
        MI computed via partial trace of the full density matrix.

        Returns:
            list of (fragment_size, average_mutual_information)
        """
        rho = self.history[time_index]
        n_env = len(self.env_indices)
        n_total = self.n_sites
        profile = []

        # System reduced density matrix (single site)
        rho_sys = rho[ix_([self.sys_idx], [self.sys_idx])]

        rng = np.random.RandomState(42)

        for f_size in range(1, n_env + 1):
            mis = []
            # Sample random subsets (up to min(n_samples, C(n_env, f_size)))
            n_combos = n_env if f_size == 1 else min(20, n_env * (n_env - 1) // 2)
            for _ in range(n_combos):
                frag = list(rng.choice(self.env_indices, f_size, replace=False))
                frag.sort()

                # Fragment reduced density matrix
                frag_mask = np.array(frag)
                rho_frag = rho[ix_(frag_mask, frag_mask)]
                S_frag = self.von_neumann_entropy(rho_frag)

                # System + fragment joint density matrix
                combined = [self.sys_idx] + frag
                rho_joint = rho[ix_(combined, combined)]
                S_joint = self.von_neumann_entropy(rho_joint)

                # MI = S(sys) + S(frag) - S(sys,frag)
                S_sys = self.von_neumann_entropy(rho_sys)
                mi = max(0.0, S_sys + S_frag - S_joint)
                mis.append(mi)

            profile.append((f_size, float(np.mean(mis))))

        return profile

    def mi_between_sites(self, i, j, time_index):
        """Quantum Mutual Information between two sites (single-exciton subspace).

        For a single-exciton system in the site basis:
        - ρ_i is 2x2: [[p_i, 0], [0, 1-p_i]]
        - ρ_j is 2x2: [[p_j, 0], [0, 1-p_j]]
        - ρ_ij is the 2x2 submatrix: [[p_i, c_ij], [c_ji, p_j]]
        where c_ij = rho[i,j] is the coherence between sites i and j.

        I(i:j) = S(ρ_i) + S(ρ_j) - S(ρ_ij)
        """
        rho = self.history[time_index]
        p_i = np.real(rho[i, i])
        p_j = np.real(rho[j, j])
        c_ij = rho[i, j]

        # Single-site reduced matrices
        rho_i = np.array([[p_i, 0], [0, 1 - p_i]])
        rho_j = np.array([[p_j, 0], [0, 1 - p_j]])

        # Two-site submatrix (single-exciton subspace)
        # Normalize by trace to make entropy comparable to single-site entropies
        rho_ij = np.array([[p_i, c_ij], [np.conj(c_ij), p_j]])
        tr_ij = np.real(np.trace(rho_ij))
        if tr_ij > 1e-15:
            rho_ij /= tr_ij

        return max(0.0, self.von_neumann_entropy(rho_i) +
                   self.von_neumann_entropy(rho_j) -
                   self.von_neumann_entropy(rho_ij))

    def redundancy_at_delta(self, delta, time_index):
        """R(delta) = number of environment fragments with I >= delta.

        R(delta) quantifies how many independent observers can each
        learn delta bits about the system.
        """
        profile = self.compute_redundancy_curve(time_index)
        return sum(1 for _, mi in profile if mi >= delta)

    def find_pointer_states(self, time_index):
        """Identify pointer states via einselection.

        Pointer states are eigenstates of the environmental monitoring
        operator — they survive decoherence and become classical.

        Computed as eigenstates of the system density matrix projected
        into the pointer basis via the environment-as-reference approach.
        """
        rho = self.history[time_index]
        env_mask = ix_(self.env_indices, self.env_indices)
        rho_env = rho[env_mask]

        # Diagonalize environment to find the einselected basis
        ev, evecs = linalg.eigh(rho_env)
        
        # The pointer states are the eigenvectors of the environmental
        # density matrix — these are the states that are robust against
        # further decoherence.
        return evecs[:, ::-1], ev[::-1]  # sorted by eigenvalue, descending

    def compute_sbs_fidelity(self, time_index):
        """Spectrum Broadcast Structure fidelity.

        SBS is achieved when I(system : fragment_k) = I(system : environment)
        for all fragments k — meaning each fragment carries complete
        classical information about the system.

        Returns SBS fidelity F in [0, 1].
        """
        rho = self.history[time_index]
        rho_sys = rho[ix_([self.sys_idx], [self.sys_idx])]
        S_sys = self.von_neumann_entropy(rho_sys)

        # I(system : all environment)
        all_env = [self.sys_idx] + self.env_indices
        rho_joint = rho[ix_(all_env, all_env)]
        rho_env = rho[ix_(np.array(self.env_indices), np.array(self.env_indices))]
        I_total = max(0.0, S_sys + self.von_neumann_entropy(rho_env) - self.von_neumann_entropy(rho_joint))

        # Average I(system : single environment site)
        I_single = 0.0
        for idx in self.env_indices:
            I_single += self.mi_between_sites(self.sys_idx, idx, time_index)
        I_single /= len(self.env_indices)

        return min(1.0, I_single / I_total) if I_total > 0 else 0.0

    def redundancy_trajectory(self, delta=0.5, step=10):
        """Time evolution of R(delta)."""
        R_t = []
        for t_idx in range(0, len(self.history), step):
            R_t.append(self.redundancy_at_delta(delta, t_idx))
        return R_t

    def print_report(self):
        """Comprehensive Darwinism report."""
        print("=" * 60)
        print("  QUANTUM DARWINISM ANALYSIS — FMO 7-site")
        print("  System: site 1  |  Environment: sites 2-7")
        print("=" * 60)

        # PIC at three time points
        for t_label, t_idx in [("Early (t=0.5 ps)", 50),
                                ("Mid (t=2.0 ps)", 200),
                                ("Late (t=10 ps)", min(999, len(self.history) - 1))]:
            if t_idx >= len(self.history):
                continue
            profile = self.compute_redundancy_curve(t_idx)
            print(f"\n  PARTIAL INFORMATION CURVE — {t_label}")
            print(f"  {'Fragment size':<16} {'I(sys:frag)':<16}")
            print(f"  {'-'*32}")
            for size, mi in profile:
                marker = " *" if mi >= profile[-1][1] * 0.5 else ""
                print(f"  {size:<16} {mi:<16.4f}{marker}")
            R_delta = self.redundancy_at_delta(0.5, t_idx)
            print(f"  R(delta=0.5): {R_delta} fragments")

        # Pointer states
        t_late = min(len(self.history) - 1, 500)
        ps, vals = self.find_pointer_states(t_late)
        print(f"\n  POINTER STATES (einselection at t=5ps)")
        print(f"  Dominant eigenvalues of rho_env:")
        for i, v in enumerate(vals[:5]):
            print(f"    lambda_{i+1} = {v:.6f}")

        # Pairwise QMI between system and each environment site
        print(f"\n  PAIRWISE QMI (site 1 with each environment site)")
        for env_idx in self.env_indices:
            mi = self.mi_between_sites(self.sys_idx, env_idx, t_late)
            print(f"    I(site 1 : site {env_idx+1}) = {mi:.4f} bits")
        sbs = self.compute_sbs_fidelity(t_late)
        print(f"\n  SPECTRUM BROADCAST STRUCTURE")
        print(f"  SBS fidelity at t=5ps: {sbs:.4f}")
        if sbs > 0.5:
            print(f"  -> Environment approaching SBS (classical limit)")
        else:
            print(f"  -> Quantum correlations persist; incomplete classicality")


if __name__ == "__main__":
    print("  Running FMO Lindblad for Darwinism analysis...")
    solver = FmoLindbladSolver()

    # Run at ENAQT-optimal dephasing (175 cm^-1)
    eff, hist = solver.run_time_evolution(dephasing_rate=175 * 0.0188, t_max=10.0)

    qd = QuantumDarwinism(hist, system_site=0)
    qd.print_report()

# API compatibility alias
FmoDarwinismAnalyzer = QuantumDarwinism
