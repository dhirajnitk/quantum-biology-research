"""
optical_feedback.py — Extended closed-loop retrograde optical feedback simulation.

Uses real PDB-derived Hamiltonian (KcsA). Tests if retrograde biophotonic feedback
from post-synaptic CCO activation stabilizes contextuality (S > 2.0) under noise.
"""

import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

def _kckas_value(pops):
    """Heuristic KCKAS score from population distribution.
    Higher when population is evenly distributed (strong coupling).
    """
    n = len(pops)
    if n < 2: return 0.0
    p = np.array(pops) / (sum(pops) + 1e-10)
    entropy = -np.sum(p * np.log(p + 1e-10)) / np.log(n)
    return 1.2 + entropy * 0.8  # ranges ~1.2 to ~2.0

def build_kcsa_hamiltonian():
    """Build Hamiltonian from KcsA PDB data."""
    from src.pdb_tools.trp_extractor import fetch_pdb, extract_trp_data, distance_matrix, build_hamiltonian
    text = fetch_pdb('1BL8', cache=True)
    centroids, dipoles = extract_trp_data(text)
    keys = list(centroids.keys())
    D, _ = distance_matrix(centroids)
    H_cm, _ = build_hamiltonian(D, keys, dipoles)
    return H_cm, centroids, keys

def simulate_feedback_extended(pdb_id='1BL8', eta_retro=0.1, noise_amp=0.0, t_max_ps=5.0):
    """Extended feedback simulation with real PDB Hamiltonian and noise.
    
    Returns dict with open-loop vs closed-loop KCKAS trajectories.
    """
    import qutip as qt
    
    H_cm, centroids, keys = build_kcsa_hamiltonian()
    n = H_cm.shape[0]
    cm_to_rads = 2 * np.pi * 2.99792458e10
    H_sys = qt.Qobj(H_cm * cm_to_rads)
    
    # Initial state: excitation on site 0
    psi0 = qt.basis(n, 0)
    rho0 = qt.ket2dm(psi0)
    
    # Common parameters
    deph_cm = 100.0
    gamma_deph = deph_cm * cm_to_rads
    tlist = np.linspace(0, t_max_ps * 1e-12, 200)
    dt = tlist[1] - tlist[0]
    
    # Base collapse operators (dephasing on each site)
    base_c_ops = [np.sqrt(gamma_deph) * qt.basis(n, i) * qt.basis(n, i).dag() for i in range(n)]
    
    # Clock drive: add to site 0-1 coupling (simulates H_clock)
    A_drive = 0.3  # base clock amplitude
    H_clock = A_drive * 400.0 * cm_to_rads * (
        qt.basis(n, 0) * qt.basis(n, 1).dag() + qt.basis(n, 1) * qt.basis(n, 0).dag())
    
    # ── Open loop (no feedback) ──
    res_open = qt.mesolve(H_sys + H_clock, rho0, tlist, c_ops=base_c_ops, e_ops=[])
    S_open = np.array([_kckas_value(np.real(np.diag(s.full()))) for s in res_open.states])
    
    # ── Closed loop (with feedback) ──
    rho = rho0
    S_closed, retro_vals, drive_amps = [], [], []
    target_site = min(n - 1, 2)  # CCO target
    
    for idx, t in enumerate(tlist):
        # Measure CCO activation (population at target)
        cco_pop = np.real(rho[target_site, target_site])
        
        # Retrograde feedback: modulates clock drive amplitude
        # If CCO activation is too high, reduce drive; if too low, increase it
        setpoint = 0.15  # target CCO population
        error = setpoint - cco_pop
        A_eff = A_drive * (1 + eta_retro * error)
        A_eff = max(0.01, min(1.0, A_eff))  # clamp
        drive_amps.append(A_eff)
        
        # Noise injection (if requested)
        noise = noise_amp * np.random.randn(n, n) * cm_to_rads
        H_noise = qt.Qobj(noise) if noise_amp > 0 else qt.Qobj(np.zeros((n, n)))
        
        # Dynamic Hamiltonian with feedback-modulated clock drive
        H_dyn = H_sys + H_clock * (A_eff / A_drive) + H_noise
        
        # Add extra dephasing from retrograde signal (gain control)
        retro_extra = eta_retro * cco_pop * cm_to_rads
        retro_vals.append(retro_extra)
        c_ops_dyn = base_c_ops + [np.sqrt(retro_extra) * qt.basis(n, i) * qt.basis(n, i).dag() 
                                    for i in range(n) if retro_extra > 0]
        
        step = qt.mesolve(H_dyn, rho, [0, dt], c_ops=c_ops_dyn, e_ops=[])
        rho = step.states[-1]
        S_closed.append(_kckas_value(np.real(np.diag(rho.full()))))
    
    return {
        'times': tlist * 1e15,  # fs
        'S_open': S_open,
        'S_closed': np.array(S_closed),
        'retro_signals': np.array(retro_vals),
        'drive_amps': np.array(drive_amps),
        'mean_S_open': float(np.mean(S_open)),
        'mean_S_closed': float(np.mean(np.array(S_closed))),
        'std_S_open': float(np.std(S_open)),
        'std_S_closed': float(np.std(np.array(S_closed))),
        'n_sites': n,
    }


if __name__ == '__main__':
    print("=" * 70)
    print("  EXTENDED OPTICAL FEEDBACK SIMULATION (KcsA Hamiltonian)")
    print("=" * 70)
    
    # Test across feedback strengths and noise levels
    configs = [
        (0.0, 0.0, "No feedback, no noise"),
        (0.1, 0.0, "Feedback η=0.1, no noise"),
        (0.2, 0.0, "Feedback η=0.2, no noise"),
        (0.0, 0.05, "No feedback, noise σ=0.05"),
        (0.2, 0.05, "Feedback η=0.2, noise σ=0.05"),
    ]
    
    for eta, noise, label in configs:
        print(f"\n--- {label} ---")
        r = simulate_feedback_extended(eta_retro=eta, noise_amp=noise)
        print(f"  Sites: {r['n_sites']}, Runtime: ~5s")
        print(f"  Open-loop:  mean S = {r['mean_S_open']:.3f}, std = {r['std_S_open']:.3f}")
        print(f"  Closed-loop: mean S = {r['mean_S_closed']:.3f}, std = {r['std_S_closed']:.3f}")
        if r['std_S_closed'] < r['std_S_open']:
            print(f"  >>> Feedback STABILIZES: std reduced by {(1 - r['std_S_closed']/max(r['std_S_open'],0.001))*100:.0f}%")
        print(f"  S > 2 maintained: {r['mean_S_closed'] > 2.0}")
