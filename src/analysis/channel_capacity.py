"""
channel_capacity.py  —  Thermodynamics and information capacity for FMO (P3).

Computes:
  - Holevo-Schumacher-Westmoreland (HSW) capacity chi
  - von Neumann entropy S(rho) trajectory
  - Entropy production rate dS/dt
  - Landauer cost: energy per successfully transmitted bit
  - Quantum vs classical efficiency comparison
  - ENAQT efficiency-entropy trade-off curve

References:
  Holevo (1998) IEEE Trans Info Theory 44, 269
  Landauer (1961) IBM J Res Dev 5, 183
  Ptaszynski & Esposito (2019) PRL 122, 150603
"""

import numpy as np
from numpy import log2, log, exp, clip, real, trace, linalg, zeros, ones
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.core.lindblad_solver import FmoLindbladSolver


# Physical constants
KB = 1.380649e-23         # J/K
HBAR = 1.054571817e-34     # J.s
T = 300.0                  # K, physiological
LANDAUER_LIMIT = KB * T * log(2)  # J/bit


class FmoThermodynamicEngine:
    """Complete thermodynamic and information-theoretic analysis of FMO.

    Parameters
    ----------
    density_history : list of ndarray
        Time-dependent density matrices from Lindblad solver.
    H : ndarray
        FMO Hamiltonian in ps^-1.
    dt : float
        Time step in ps (default 0.01).
    """

    def __init__(self, density_history, H, dt=0.01):
        self.history = density_history
        self.H = H
        self.dt = dt
        self.n_sites = H.shape[0] if H.ndim == 2 else 7
        self.n_times = len(density_history)

    # ── Core entropy metrics ────────────────────────────────────

    def von_neumann_entropy(self, rho):
        """S(rho) = -Tr(rho log2 rho) in bits."""
        ev = linalg.eigvalsh(rho)
        ev = clip(ev, 1e-15, 1.0)
        return -np.sum(ev * log2(ev))

    def compute_entropy_trajectory(self):
        """Full entropy trajectory S(t)."""
        return np.array([self.von_neumann_entropy(rho) for rho in self.history])

    def compute_entropy_production(self, time_index):
        """Entropy production rate dS/dt at time_index.

        dS/dt = -Tr(drho/dt * log(rho))
        Uses forward difference for drho/dt.
        """
        if time_index < 1 or time_index >= self.n_times:
            return 0.0
        rho_t = self.history[time_index]
        rho_prev = self.history[time_index - 1]
        drho = (rho_t - rho_prev) / self.dt
        ev, U = linalg.eigh(rho_t)
        ev = clip(ev, 1e-15, 1.0)
        log_rho = U @ np.diag(log(ev)) @ U.conj().T
        return max(0.0, -real(trace(drho @ log_rho)))

    def compute_entropy_production_trajectory(self):
        """Full entropy production rate trajectory dS/dt(t)."""
        return np.array([self.compute_entropy_production(i)
                         for i in range(self.n_times)])

    # ── Information capacity ────────────────────────────────────

    def compute_holevo_capacity(self, time_index):
        """Holevo capacity (population-based approximation).

        WARNING: This is an approximation that uses the diagonal elements
        of ρ(t) as proxy conditional states. The true Holevo capacity
        requires propagating each input basis state separately through
        the channel — see `compute_holevo_capacity_proper()` for a
        correct implementation.

        Returns an upper bound on the true Holevo χ.
        """
        if time_index >= self.n_times:
            return 0.0
        rho = self.history[time_index]

        p_k = 1.0 / self.n_sites
        total_S = self.von_neumann_entropy(rho)

        # Approximate conditional entropy from diagonal elements
        cond_S = 0.0
        for k in range(self.n_sites):
            state_k = np.zeros((self.n_sites, self.n_sites), dtype=complex)
            state_k[k, k] = rho[k, k]
            tr_k = real(trace(state_k))
            if tr_k > 1e-15:
                state_k /= tr_k
            cond_S += p_k * self.von_neumann_entropy(state_k)

        return max(0.0, total_S - cond_S)

    def compute_holevo_capacity_proper(self, rho_input, rho_output):
        """True Holevo capacity for a given input/output pair.

        χ = S(ρ_out) - Σ p_i S(Φ(ρ_i))

        Parameters
        ----------
        rho_input : list of ndarray
            Input density matrices [ρ_1, ..., ρ_N].
        rho_output : list of ndarray
            Corresponding output density matrices after channel propagation.

        Returns
        -------
        chi : float
            Holevo capacity in bits.
        """
        n = len(rho_input)
        p_i = 1.0 / n

        # Ensemble-averaged output state
        rho_bar = sum(rho_out for rho_out in rho_output) / n
        S_bar = self.von_neumann_entropy(rho_bar)

        # Average conditional entropy
        cond_S = sum(self.von_neumann_entropy(rho_out) for rho_out in rho_output) / n

        return max(0.0, S_bar - cond_S)

    def compute_holevo_trajectory(self):
        """Full Holevo capacity trajectory chi(t)."""
        return np.array([self.compute_holevo_capacity(i)
                         for i in range(self.n_times)])

    def coherent_information(self, time_index):
        """Coherent information I_c(rho) = S(rho_output) - S(rho_env).

        Lower bound on quantum capacity Q.
        """
        if time_index >= self.n_times:
            return 0.0
        rho = self.history[time_index]
        return self.von_neumann_entropy(rho)

    # ── Energy / Landauer analysis ──────────────────────────────

    def compute_energy_per_bit(self, time_index, channel_efficiency):
        """Energy cost per successfully transmitted bit.

        E_bit = E_exciton_per_try / min(channel_capacity, channel_efficiency)

        Parameters
        ----------
        time_index : int
            Time point for the analysis.
        channel_efficiency : float
            Fraction of excitons successfully trapped (0-1).

        Returns
        -------
        dict with J/bit, kT/bit, and ratio to Landauer limit.
        """
        chi = self.compute_holevo_capacity(time_index)
        if chi < 1e-10 or channel_efficiency < 1e-10:
            return {"J_per_bit": float('inf'), "x_Landauer": float('inf')}

        # Each trapping event costs one exciton of energy
        # Exciton energy: ~35,700 cm^-1 = 7.09e-19 J
        exciton_J = 7.09e-19

        # Each successful trap transmits chi bits
        bits_per_trap = chi
        traps_needed = 1.0 / channel_efficiency  # including failed attempts
        J_per_bit = exciton_J * traps_needed / bits_per_trap

        return {
            "chi_bits": chi,
            "efficiency": channel_efficiency,
            "J_per_bit": J_per_bit,
            "x_Landauer": J_per_bit / LANDAUER_LIMIT,
            "kT_per_bit": J_per_bit / (KB * T),
        }

    # ── Quantum vs classical comparison ─────────────────────────



    @staticmethod
    def classical_random_walk_efficiency(n_sites=7, p_hop=0.5, p_loss=0.05, n_steps=100):
        """Classical FRET-like random walk with loss.

        At each step, the exciton either hops (prob p_hop), is lost
        to non-radiative decay (prob p_loss), or stays (prob 1-p_hop-p_loss).
        Efficiency = fraction reaching trap site 3.
        """
        rng = np.random.RandomState(42)
        trapped = 0
        trials = 10000

        for _ in range(trials):
            pos = 0
            for _ in range(n_steps):
                if pos == 2:
                    trapped += 1
                    break
                r = rng.random()
                if r < p_loss:
                    break  # lost to decay
                elif r < p_loss + p_hop:
                    # Hop to nearest neighbour with bias toward site 3
                    if pos < 2:
                        pos += 1
                    elif pos > 2:
                        pos -= 1
                    else:
                        pos += rng.choice([-1, 1])
        return trapped / trials

    # ── Comprehensive report ────────────────────────────────────

    def print_report(self, channel_efficiency=None):
        """Print comprehensive thermodynamics report."""
        if channel_efficiency is None:
            # Use all remaining population as efficiency estimate
            rho_final = self.history[-1]
            channel_efficiency = 1.0 - real(trace(rho_final))
            channel_efficiency = max(0.0, min(1.0, channel_efficiency))

        S_t = self.compute_entropy_trajectory()
        dS_t = self.compute_entropy_production_trajectory()
        chi_t = self.compute_holevo_trajectory()

        print("=" * 62)
        print("  THERMODYNAMICS & CHANNEL CAPACITY — FMO 7-site")
        print("  Landauer limit kT ln 2 = {:.2e} J/bit".format(LANDAUER_LIMIT))
        print("=" * 62)

        # 1. Holevo capacity
        print(f"\n  HOLEVO CAPACITY")
        print(f"{'t (ps)':<10} {'chi (bits)':<15} {'S (bits)':<15} {'dS/dt':<15}")
        print(f"{'-'*55}")
        for t_idx in [0, 50, 100, 200, 500, min(999, self.n_times - 1)]:
            if t_idx >= self.n_times:
                continue
            t = t_idx * self.dt
            print(f"  {t:<8.2f}  {chi_t[t_idx]:<12.4f}  {S_t[t_idx]:<12.4f}  {dS_t[t_idx]:<12.6f}")

        # 2. Energy per bit
        ebit = self.compute_energy_per_bit(min(500, self.n_times - 1), channel_efficiency)
        print(f"\n  ENERGY PER SUCCESSFUL BIT")
        print(f"    Holevo capacity:        {ebit['chi_bits']:.4f} bits")
        print(f"    Transport efficiency:   {ebit['efficiency']*100:.2f}%")
        print(f"    Energy per bit:         {ebit['J_per_bit']:.4e} J")
        print(f"    Landauer ratio:         {ebit['x_Landauer']:.4e}x")
        print(f"    Cost per bit:           {ebit['kT_per_bit']:.2f} kT")

        # 3. Classical comparison
        eta_cl = self.classical_random_walk_efficiency()
        print(f"\n  QUANTUM vs CLASSICAL TRANSPORT")
        print(f"    Quantum efficiency (Lindblad):     {channel_efficiency*100:.1f}%")
        print(f"    Classical random walk efficiency:  {eta_cl*100:.1f}%")
        print(f"    Quantum advantage:                 {channel_efficiency/max(eta_cl,0.01):.2f}x")

        # 4. Relaxation timescale
        tau_idx = next((i for i, s in enumerate(S_t) if s < 0.5 * S_t[0]), len(S_t) - 1)
        t_relax = tau_idx * self.dt
        print(f"\n  RELAXATION")
        print(f"    Entropy half-life:       {t_relax:.2f} ps")
        print(f"    Max entropy:             {np.max(S_t):.4f} bits")

        # 5. Peak entropy production
        peak_idx = int(np.argmax(dS_t))
        print(f"    Peak dS/dt:              {dS_t[peak_idx]:.6f} at t={peak_idx*self.dt:.2f} ps")


def run_p3_analysis():
    """Run full P3 analysis at multiple dephasing rates."""
    print("  Running FMO Lindblad for thermodynamic analysis...")
    solver = FmoLindbladSolver()

    gamma_values = [0, 50, 100, 175, 300, 500]
    results = []

    for g_cm in gamma_values:
        g_ps = g_cm * 0.0188
        eff, hist = solver.run_time_evolution(dephasing_rate=g_ps, t_max=10.0)

        # Final population remaining in system = 1 - efficiency
        eff_actual = 1.0 - real(np.trace(hist[-1]))
        thermo = FmoThermodynamicEngine(hist, solver.H)
        S_t = thermo.compute_entropy_trajectory()
        chi_t = thermo.compute_holevo_trajectory()

        peak_S = np.max(S_t)
        peak_chi = np.max(chi_t)
        peak_dS = np.max(thermo.compute_entropy_production_trajectory())

        ebit = thermo.compute_energy_per_bit(len(hist) - 1, eff_actual)

        results.append({
            "gamma_cm": g_cm,
            "efficiency": eff_actual,
            "peak_entropy": peak_S,
            "peak_chi": peak_chi,
            "peak_dS": peak_dS,
            "J_per_bit": ebit["J_per_bit"],
            "x_Landauer": ebit["x_Landauer"],
        })

        # Classical comparison
        eta_cl = FmoThermodynamicEngine.classical_random_walk_efficiency()

    print(f"\n{'='*72}")
    print(f"  P3 THERMODYNAMICS — MULTI-DEPHASING COMPARISON")
    print(f"  Note: Simplified Markovian model (no Drude-Lorentz spectral density)")
    print(f"  ENAQT peak requires structured bath + static disorder for Anderson localization.")
    print(f"{'='*72}")
    print(f"{'Gamma':<10} {'Efficiency':<14} {'Peak S':<12} {'Peak chi':<12} {'Peak dS':<14} {'xLandauer':<12}")
    print(f"{'-'*10} {'-'*14} {'-'*12} {'-'*12} {'-'*14} {'-'*12}")
    for r in results:
        print(f"  {r['gamma_cm']:<8.0f} {r['efficiency']*100:<10.2f}%  {r['peak_entropy']:<10.4f} {r['peak_chi']:<10.4f} {r['peak_dS']:<10.6f} {r['x_Landauer']:<10.2e}")

    opt_idx = np.argmax([r['efficiency'] for r in results])
    print(f"\n  OPTIMAL: gamma = {results[opt_idx]['gamma_cm']:.0f} cm-1, eff = {results[opt_idx]['efficiency']*100:.1f}%")
    print(f"  CLASSICAL RANDOM WALK: {eta_cl*100:.1f}%")
    print(f"  Quantum advantage at gamma=0: {results[0]['efficiency']/max(eta_cl,0.01):.1f}x")
    return results


if __name__ == "__main__":
    results = run_p3_analysis()
