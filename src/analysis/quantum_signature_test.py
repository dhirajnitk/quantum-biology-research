"""
quantum_signature_test.py  —  Experimental tests to prove quantum origin
of the Z-channel gating mechanism.

Tests:
  1. Superradiance scaling: classical (rate ∝ N) vs quantum (rate ∝ N²)
  2. Phase-modulation interference: classical (phase-independent) vs
     quantum (cos²θ dependence)
  3. Bell-type inequality test: theoretical framework for biological systems

If the gating follows classical statistics, the spatial ensemble model
predicts linear scaling. If quantum superradiance is at work, the
per-core emission rate itself increases with N (Dicke superradiance),
producing a characteristic quadratic-to-linear crossover.

References
----------
[Dicke1954] Dicke (1954). "Coherence in spontaneous radiation processes."
[Babcock2024] Babcock et al. (2024) JPCB — Trp superradiance confirmed.
[Bell1964] Bell (1964). "On the Einstein Podolsky Rosen paradox."
"""

import numpy as np
from numpy import exp, sqrt, cos, sin, pi, log2, random


# ═══════════════════════════════════════════════════════════════════
#  1. Superradiance scaling test
# ═══════════════════════════════════════════════════════════════════

def classical_emission_rate(n_emitters, p_single=7.79e-4):
    """Classical: emitters are independent, rate scales linearly.

    Each emitter fires independently with probability p_single.
    Total expected photons = N * p_single.
    """
    return n_emitters * p_single


def superradiant_emission_rate(n_emitters, p_single=7.79e-4,
                                cooperative_factor=1.0):
    """Quantum superradiant: emitters are phase-coherent.

    In the Dicke superradiant regime, the emission rate scales
    quadratically for small N: rate ∝ N².
    The cooperative factor depends on the geometric arrangement
    and dipole-dipole coupling strength.

    For a fully coherent ensemble:
        rate = N * p_single * (1 + (N-1) * μ)

    where μ is the cooperative enhancement factor (0 ≤ μ ≤ 1).
    μ = 0 → classical independent emitters
    μ = 1 → fully coherent superradiance (rate ∝ N² for small N)

    [Dicke1954] Eq. 4.1, [Babcock2024] Fig 3.
    """
    rate_per_emitter = p_single * (1 + (n_emitters - 1) * cooperative_factor)
    return n_emitters * max(rate_per_emitter, 0)


def superradiance_crossover(n_max=1000, cooperative_factor=1.0):
    """Compute classical vs quantum predictions and the ratio.

    The ratio quantum/classical is the key observable:
    - ratio = 1 → classical
    - ratio > 1 → superradiant enhancement
    - ratio ∝ N for small N in fully coherent regime
    """
    ns = np.arange(1, n_max + 1)
    classical = np.array([classical_emission_rate(n) for n in ns])
    quantum = np.array([superradiant_emission_rate(n, cooperative_factor=cooperative_factor)
                        for n in ns])
    ratio = quantum / classical
    return ns, classical, quantum, ratio


def detect_superradiance_threshold(ns, ratio, threshold=1.5):
    """Find the ensemble size at which quantum enhancement becomes detectable.

    Returns the smallest N where ratio > threshold.
    """
    idx = np.where(ratio > threshold)[0]
    if len(idx) > 0:
        return ns[idx[0]]
    return None


# ═══════════════════════════════════════════════════════════════════
#  2. Phase-modulation interference test
# ═══════════════════════════════════════════════════════════════════

def classical_phase_response(theta, amplitude=1.0):
    """Classical prediction: phase modulation has no effect on intensity.

    Incoherent classical emitters respond only to intensity, not phase.
    """
    return amplitude * np.ones_like(theta)


def quantum_phase_response(theta, amplitude=1.0, visibility=1.0):
    """Quantum prediction: interference fringe from phase-coherent emitters.

    For a quantum coherent ensemble, the emission probability follows:
        I(θ) = amplitude * [1 + V * cos(θ)]

    where V is the visibility (0 ≤ V ≤ 1).
    V = 0 → classical (no interference)
    V = 1 → fully coherent (maximum interference)

    This is the signature of quantum phase coherence.
    """
    return amplitude * (1.0 + visibility * cos(theta))


def fringe_visibility(classical_response, quantum_response):
    """Compute fringe visibility V = (I_max - I_min) / (I_max + I_min).

    V = 0 → classical (no fringes)
    V > 0 → quantum coherence detected
    """
    I_min = np.min(quantum_response)
    I_max = np.max(quantum_response)
    return (I_max - I_min) / (I_max + I_min)


# ═══════════════════════════════════════════════════════════════════
#  3. Bell-type inequality for biological systems
# ═══════════════════════════════════════════════════════════════════

def chsh_inequality_classical_bound():
    """Clauser-Horne-Shimony-Holt (CHSH) inequality.

    For any local hidden variable theory:
        |S| ≤ 2

    where S = E(a,b) + E(a,b') + E(a',b) - E(a',b')
    and E is the correlation function.

    Quantum mechanics predicts |S| ≤ 2√2 ≈ 2.828.
    Violation of |S| ≤ 2 proves non-classical correlations.
    """
    return 2.0


def chsh_quantum_bound():
    """Tsirelson's bound: maximum quantum violation of CHSH."""
    return 2.0 * sqrt(2.0)


def simulate_chsh_test(n_measurements=10000, entanglement_fidelity=0.5):
    """Simulate a CHSH test on two spatially separated Trp ensembles.

    Parameters
    ----------
    n_measurements : int
        Number of coincidence measurements.
    entanglement_fidelity : float
        How close the two ensembles are to a maximally entangled state.
        0 = completely mixed (classical), 1 = maximally entangled.

    Returns
    -------
    S : float
        CHSH correlation value. S > 2 indicates non-classical.
    """
    # Measurement angles (standard CHSH settings)
    angles_a = [0, pi / 2]
    angles_b = [pi / 4, 3 * pi / 4]

    correlations = np.zeros((2, 2))
    for i, a in enumerate(angles_a):
        for j, b in enumerate(angles_b):
            # Quantum correlation for a partially entangled state
            # E(a,b) = entanglement_fidelity * cos(2*(a-b))
            # For max entanglement: E = cos(2*(a-b))
            # Classical limit: E = entanglement_fidelity * something (≤ 2)
            e = entanglement_fidelity * cos(2 * (a - b))
            correlations[i, j] = e

    S = correlations[0, 0] + correlations[0, 1] + correlations[1, 0] - correlations[1, 1]
    return S


# ═══════════════════════════════════════════════════════════════════
#  4. Experimental design
# ═══════════════════════════════════════════════════════════════════

class QuantumSignatureExperiment:
    """Design specification for an experiment to prove quantum origin."""

    @staticmethod
    def test_1_superradiance():
        """Test 1: Emission rate scaling.

        Stimulate membrane patches with increasing UV laser intensity
        (N emitters). Measure the fluorescence intensity I(N).

        Classical: I ∝ N
        Quantum:   I ∝ N² for small N, crossing over to I ∝ N at large N

        Observable: ratio I(N) / (N * I_single) > 1 for small N
        """
        return {
            "name": "Superradiance scaling test",
            "stimulus": "UV laser pulses (327 nm), increasing power",
            "measurement": "Time-resolved fluorescence intensity",
            "classical_prediction": "I(N) / (N * I_1) = 1 (flat)",
            "quantum_prediction": "I(N) / (N * I_1) > 1 for small N, decaying to 1 at large N",
            "falsifiable": "If ratio = 1 for all N, system is classical.",
        }

    @staticmethod
    def test_2_phase_interference():
        """Test 2: Phase-modulation interference.

        Split a UV laser into two paths with variable phase delay θ.
        Recombine on the membrane patch. Measure gating current.

        Classical: I(θ) = constant (phase-independent)
        Quantum:   I(θ) ∝ cos²(θ) (interference fringes)

        Observable: fringe visibility V = (I_max - I_min) / (I_max + I_min)
        """
        return {
            "name": "Phase interference test",
            "stimulus": "Phase-modulated UV laser (two-path interferometer)",
            "measurement": "Membrane ionic current vs phase angle",
            "classical_prediction": "I(θ) = constant (V = 0)",
            "quantum_prediction": "I(θ) ∝ 1 + V cos(θ) with V > 0",
            "falsifiable": "If V = 0 within measurement precision, system is classical.",
        }

    @staticmethod
    def test_3_chsh_inequality():
        """Test 3: CHSH-Bell inequality.

        Requires two spatially separated Trp ensembles in the same cell,
      each driven by a phase-coherent UV laser at independently
        controllable angles a, b.

        Measure coincidence gating events: both gates open simultaneously.
        Compute S = E(a,b) + E(a,b') + E(a',b) - E(a',b').

        Classical: |S| ≤ 2
        Quantum:   |S| ≤ 2√2 ≈ 2.828

        Observable: S > 2
        """
        return {
            "name": "CHSH-Bell inequality test",
            "stimulus": "Two phase-coherent UV lasers with variable angles",
            "measurement": "Coincidence gating events at two spatially separated synaptic sites",
            "classical_prediction": "S ≤ 2",
            "quantum_prediction": "S can exceed 2, up to 2√2",
            "falsifiable": "If S ≤ 2, local hidden variables cannot be ruled out.",
        }


# ═══════════════════════════════════════════════════════════════════
#  Demo
# ═══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 72)
    print("  QUANTUM SIGNATURE TESTS FOR Z-CHANNEL GATING")
    print("  Distinguishing quantum superradiance from classical noise")
    print("=" * 72)

    # ── Test 1: Superradiance scaling ──────────────────────────
    print(f"\n  TEST 1: SUPERRADIANCE SCALING")
    ns, classical, quantum, ratio = superradiance_crossover(200, cooperative_factor=0.5)
    print(f"  {'N':<8} {'Classical':<14} {'Quantum':<14} {'Ratio':<10}")
    print(f"  {'-'*46}")
    for i in [0, 9, 49, 99, 199]:
        print(f"  {ns[i]:<8} {classical[i]:<14.4f} {quantum[i]:<14.4f} {ratio[i]:<10.4f}")

    threshold_n = detect_superradiance_threshold(ns, ratio, 1.5)
    if threshold_n:
        print(f"\n  -> Superradiance detectable at N > {threshold_n} (ratio > 1.5x)")
    print(f"  -> Falsifiable: measure I(N)/I(1) ratio. Classical flat at 1.")

    # ── Test 2: Phase interference ─────────────────────────────
    print(f"\n  TEST 2: PHASE INTERFERENCE")
    theta = np.linspace(0, 2 * pi, 5)
    classical_i = classical_phase_response(theta)
    quantum_i = quantum_phase_response(theta, visibility=0.8)
    print(f"  {'theta':<10} {'Classical':<14} {'Quantum':<14}")
    print(f"  {'-'*38}")
    for i in range(len(theta)):
        print(f"  {theta[i]:<10.2f} {classical_i[i]:<14.4f} {quantum_i[i]:<14.4f}")
    v = fringe_visibility(classical_i, quantum_i)
    print(f"\n  -> Fringe visibility: V = {v:.4f} (V=0 is classical)")

    # ── Test 3: CHSH inequality ────────────────────────────────
    print(f"\n  TEST 3: CHSH-BELL INEQUALITY")
    for fidelity in [0.0, 0.3, 0.5, 0.7, 1.0]:
        S = simulate_chsh_test(entanglement_fidelity=fidelity)
        status = "QUANTUM" if S > 2 else "CLASSICAL"
        print(f"  Entanglement fidelity: {fidelity:.1f}  S = {S:.4f}  [{status}]")
    print(f"\n  Classical bound: |S| <= 2")
    print(f"  Quantum bound:   |S| <= 2*sqrt(2) ~ 2.828")

    # ── Summary ────────────────────────────────────────
    print(f"\n{'='*72}")
    print(f"  EXPERIMENTAL ROADMAP")
    print(f"{'='*72}")
    for test in [QuantumSignatureExperiment.test_1_superradiance(),
                 QuantumSignatureExperiment.test_2_phase_interference(),
                 QuantumSignatureExperiment.test_3_chsh_inequality()]:
        print(f"\n  {test['name']}")
        print(f"    Stimulus:  {test['stimulus']}")
        print(f"    Measure:   {test['measurement']}")
        print(f"    Quantum:   {test['quantum_prediction']}")
        print(f"    Classical: {test['classical_prediction']}")
        print(f"    Verdict:   {test['falsifiable']}")
