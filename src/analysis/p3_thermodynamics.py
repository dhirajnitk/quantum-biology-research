"""
p3_thermodynamics.py  —  Comprehensive thermodynamic analysis for P3 paper.

Computes across multiple dephasing rates:
  1. ENAQT efficiency curve with Anderson localization (via static disorder)
  2. von Neumann entropy trajectory S(t)
  3. Entropy production rate dS/dt
  4. Holevo capacity chi(t)
  5. Landauer energy cost per bit
  6. Quantum vs classical random walk comparison
  7. Efficiency-entropy trade-off

This completes P3 requirements for the thermodynamic cost paper.

References:
  Landauer (1961) IBM J Res Dev 5, 183
  Ptaszynski & Esposito (2019) PRL 122, 150603
  Dorfman et al. (2013) PNAS 110, 2746
"""

import numpy as np
from numpy import log2, log, real, trace, linalg, clip
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from core.lindblad_solver import FmoLindbladSolver
from analysis.channel_capacity import FmoThermodynamicEngine


KB = 1.380649e-23
T = 300.0
LANDAUER = KB * T * log(2)
EXCITON_J = 7.09e-19  # 35,700 cm^-1 in Joules


def enaqt_with_disorder(gamma_range=None, disorder_std=50.0):
    """ENAQT efficiency curve with static disorder for Anderson localization.

    With disorder: efficiency is LOW at gamma=0 (Anderson localization),
    RISES at intermediate gamma (ENAQT), DROPS at high gamma (Zeno).
    """
    if gamma_range is None:
        gamma_range = np.logspace(-1, 3, 30)  # 0.1 to 1000 cm^-1

    effs = []
    for g_cm in gamma_range:
        g_ps = g_cm * 0.0188
        solver = FmoLindbladSolver(disorder_std=disorder_std)
        eff, _ = solver.run_time_evolution(dephasing_rate=g_ps, t_max=10.0)
        effs.append(eff)

    return gamma_range, np.array(effs)


def entropy_trajectory(solver, gamma_cm, t_max=10.0):
    """Full entropy and Holevo trajectories at given dephasing."""
    g_ps = gamma_cm * 0.0188
    eff, hist = solver.run_time_evolution(dephasing_rate=g_ps, t_max=t_max)

    thermo = FmoThermodynamicEngine(hist, solver.H)
    S_t = thermo.compute_entropy_trajectory()
    dS_t = thermo.compute_entropy_production_trajectory()
    chi_t = thermo.compute_holevo_trajectory()

    return eff, np.array(S_t), np.array(dS_t), np.array(chi_t)


def run_full_p3(disorder_std=50.0):
    """Run complete P3 analysis."""
    print("=" * 68)
    print("  P3 THERMODYNAMICS — Complete Analysis")
    print(f"  Static disorder: {disorder_std} cm-1 (Anderson localization)")
    print(f"  Landauer limit: {LANDAUER:.2e} J/bit at {T} K")
    print("=" * 68)

    # 1. ENAQT curve with disorder
    print("\n  [1] ENAQT EFFICIENCY CURVE (with disorder)")
    gammas, effs = enaqt_with_disorder(disorder_std=disorder_std)
    print(f"  {'Gamma (cm-1)':<18} {'Efficiency':<15}")
    print(f"  {'-'*33}")
    for i in range(0, len(gammas), 3):
        print(f"  {gammas[i]:<16.2f} {effs[i]*100:<10.4f}%")

    opt_idx = int(np.argmax(effs))
    print(f"\n  ENAQT OPTIMAL: gamma = {gammas[opt_idx]:.1f} cm-1, eff = {effs[opt_idx]*100:.2f}%")
    print(f"  Gamma=0 efficiency:     {effs[0]*100:.2f}%")
    print(f"  Enhancement from gamma=0 to optimal: {effs[opt_idx]/max(effs[0],0.001):.2f}x")

    # 2. Detailed analysis at optimal gamma
    print(f"\n  [2] DETAILED THERMODYNAMICS AT OPTIMAL GAMMA")
    solver = FmoLindbladSolver(disorder_std=disorder_std)
    eff_gamma_opt, S_t, dS_t, chi_t = entropy_trajectory(
        solver, gammas[opt_idx])

    dt = 0.01
    print(f"  {'t (ps)':<10} {'S(t) (bits)':<16} {'dS/dt (kB/ps)':<18} {'chi(t) (bits)':<18}")
    print(f"  {'-'*62}")
    for t_idx in [0, 50, 100, 200, 500, 999]:
        if t_idx >= len(S_t):
            continue
        t = t_idx * dt
        print(f"  {t:<8.2f}  {S_t[t_idx]:<14.4f}  {dS_t[t_idx]:<16.6f}  {chi_t[t_idx]:<14.4f}")

    # 3. Landauer cost
    rho_final = None
    _, hist = solver.run_time_evolution(gammas[opt_idx] * 0.0188, t_max=10.0)
    rho_final = hist[-1]
    efficiency = 1.0 - real(np.trace(rho_final))
    peak_chi = np.max(chi_t)

    J_per_bit = EXCITON_J / max(efficiency * peak_chi, 1e-10)
    landauer_ratio = J_per_bit / LANDAUER

    print(f"\n  [3] LANDAUER ENERGY COST")
    print(f"    Peak Holevo capacity:          {peak_chi:.4f} bits")
    print(f"    Transport efficiency:           {efficiency*100:.2f}%")
    print(f"    Energy per successful bit:      {J_per_bit:.4e} J")
    print(f"    Landauer ratio:                 {landauer_ratio:.4e}x")
    print(f"    Cost per bit:                   {J_per_bit/(KB*T):.2f} kT")

    # 4. Comparison across gamma values
    print(f"\n  [4] MULTI-DEPHASING COMPARISON")
    print(f"  {'Gamma':<10} {'Efficiency':<14} {'Peak S':<12} {'Peak chi':<12} {'Peak dS':<14} {'xLandauer':<14}")
    print(f"  {'-'*10} {'-'*14} {'-'*12} {'-'*12} {'-'*14} {'-'*14}")
    for g_cm in [0, 10, 50, 100, 175, 300, 500]:
        g_ps = g_cm * 0.0188
        s = FmoLindbladSolver(disorder_std=disorder_std)
        e, hist2 = s.run_time_evolution(g_ps, t_max=10.0)
        tr = FmoThermodynamicEngine(hist2, s.H)
        S_max = np.max(tr.compute_entropy_trajectory())
        dS_max = np.max(tr.compute_entropy_production_trajectory())
        chi_max = np.max(tr.compute_holevo_trajectory())
        eff = 1.0 - real(np.trace(hist2[-1]))
        lr = EXCITON_J / max(eff * chi_max, 1e-10) / LANDAUER
        print(f"  {g_cm:<8.0f}  {eff*100:<10.2f}%  {S_max:<10.4f}  {chi_max:<10.4f}  {dS_max:<12.6f}  {lr:<12.2e}")

    # 5. Classical comparison
    print(f"\n  [5] QUANTUM vs CLASSICAL")
    cl_eff = FmoThermodynamicEngine.classical_random_walk_efficiency()
    print(f"    Quantum optimal efficiency:     {effs[opt_idx]*100:.2f}%")
    print(f"    Classical random walk:          {cl_eff*100:.1f}%")
    print(f"    Quantum enhancement:            {effs[opt_idx]/max(cl_eff,0.01):.2f}x")

    print(f"\n{'='*68}")
    print(f"  P3 SUMMARY")
    print(f"  Optimal dephasing:    {gammas[opt_idx]:.1f} cm-1")
    print(f"  Peak efficiency:      {effs[opt_idx]*100:.2f}%")
    print(f"  Landauer ratio:       {landauer_ratio:.2e}x")
    print(f"  Peak entropy prod.:   {np.max(dS_t):.6f} kB/ps")
    print(f"  Max channel capacity: {peak_chi:.4f} bits")
    print(f"{'='*68}")

    return {
        "gammas": gammas,
        "effs": effs,
        "opt_gamma": gammas[opt_idx],
        "opt_eff": effs[opt_idx],
        "landauer_ratio": landauer_ratio,
        "peak_chi": peak_chi,
        "peak_dS": np.max(dS_t),
    }


if __name__ == "__main__":
    run_full_p3(disorder_std=50.0)
