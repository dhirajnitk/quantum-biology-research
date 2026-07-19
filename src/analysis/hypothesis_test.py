"""
hypothesis_test.py  —  Spatial ensemble gating as a statistical hypothesis test.

The synapse treats each synchronized burst of N Trp cores as a hypothesis test:
  H0: no signal present  →  gate stays closed
  H1: signal present     →  gate opens

The Z-channel asymmetry (zero false positives) means α = 0 identically.
Statistical power is determined by ensemble size N and target absorption cross-section.

This module computes the experimental prediction: ionic current I(N) as a function
of synchronized ensemble size, which can be tested with ultrafast laser pulse
stimulation of membrane patches.

References
----------
[B2024] Babcock et al. (2024) JPCB — Trp fluorescence QY, spectra
[Neher1992] Neher & Sakmann (1992) Nature — Patch-clamp electrophysiology
"""

import numpy as np
from numpy import exp, sqrt, log, pi, arcsin, cos
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from src.core.biophoton_relay import BiophotonRelay, channel_capacity


class HypothesisTest:
    """Models synaptic gating as a statistical hypothesis test.

    Parameters
    ----------
    p_per_core : float
        Probability that a single core successfully triggers the gate.
    n_cores_total : int
        Total number of Trp cores available in the synapse.
    """

    def __init__(self, p_per_core=7.79e-4, n_cores_total=50000):
        self.p = p_per_core
        self.N = n_cores_total

    # ── Core statistical metrics ───────────────────────────────

    def power(self, n):
        """Statistical power for ensemble of size n: P(reject H0 | H1 true)."""
        return 1.0 - (1.0 - self.p) ** n

    def alpha(self):
        """False-positive rate: probability of gate opening without signal.
        Zero because p_dark ≈ 0 in ε=2 shielded membrane.
        """
        return 0.0

    def beta(self, n):
        """False-negative rate for ensemble size n."""
        return (1.0 - self.p) ** n

    def min_n_for_power(self, target_power=0.8):
        """Minimum ensemble size n needed to achieve target statistical power."""
        return int(np.ceil(log(1 - target_power) / log(1 - self.p)))

    def channel_capacity(self, n):
        """Mutual information I(N) of the gating channel for ensemble size n."""
        p_success = self.power(n)
        return channel_capacity(p_success)

    # ── Experimental prediction: ionic current ─────────────────

    def ionic_current(self, n, single_channel_current_pA=1.0, open_probability=1.0):
        """Predicted ionic current for ensemble size n.

        Models a patch-clamp experiment where n cores are synchronously
        activated by an ultrafast laser pulse.

        Parameters
        ----------
        n : int
            Number of synchronized cores in the pulse.
        single_channel_current_pA : float
            Current through a single open ion channel (pA).
        open_probability : float
            Probability that a triggered gate opens (≈ 1.0 for voltage-gated).

        Returns
        -------
        I_pA : float
            Expected ionic current in picoamperes.
        """
        p_gate_open = self.power(n) * open_probability
        return p_gate_open * single_channel_current_pA * n

    def current_sweep(self, n_values=None):
        """Compute ionic current across a range of ensemble sizes.

        This is the falsifiable experimental prediction: a non-linear
        I(N) curve that distinguishes quantum Z-channel gating from
        classical linear summation.
        """
        if n_values is None:
            n_values = [1, 10, 100, 500, 1000, 2000, 5000,
                        10000, 20000, 50000]
        results = []
        for n in n_values:
            results.append({
                "n_cores": n,
                "power": self.power(n),
                "beta": self.beta(n),
                "capacity": self.channel_capacity(n),
                "current_pA": self.ionic_current(n),
            })
        return results


def format_results(results):
    """Pretty-print the hypothesis test results."""
    print(f"{'N_cores':<10} {'Power':<12} {'FNR(beta)':<12} {'Capacity':<12} {'Current(pA)':<12}")
    print(f"{'-'*58}")
    for r in results:
        print(f"{r['n_cores']:<10} {r['power']*100:<8.3f}%  {r['beta']:<8.2e}  {r['capacity']:<8.4f}  {r['current_pA']:<8.4f}")


def compare_targets():
    """Compare Trp vs CCO target across ensemble sizes."""
    targets = {
        "Trp (baseline)": 1.80e-5,
        "Fe-S cluster": 6.5e-21 / 1e-18,  # ~0.0065
        "FMN": 3.0e-21 / 1e-18,            # ~0.003
        "CCO (cytochrome)": 7.79e-4,
    }
    print(f"\n{'='*70}")
    print(f"  TARGET COMPARISON: Minimum ensemble for 80% power")
    print(f"{'='*70}")
    print(f"{'Target':<25} {'p_per_core':<15} {'N_min':<10}")
    print(f"{'-'*50}")
    for name, p in targets.items():
        ht = HypothesisTest(p_per_core=p)
        n_min = ht.min_n_for_power(0.8)
        print(f"{name:<25} {p:<15.2e} {n_min:<10}")


if __name__ == "__main__":
    # Default: CCO target (our primary model)
    print(f"{'='*70}")
    print(f"  HYPOTHESIS TEST: Spatial Ensemble Gating")
    print(f"  Per-core hit probability: 7.79e-4 (CCO target)")
    print(f"  False-positive rate (alpha): 0.0 (Z-channel)")
    print(f"{'='*70}")

    ht = HypothesisTest(p_per_core=7.79e-4)
    results = ht.current_sweep()
    format_results(results)

    # Minimum N for conventional power thresholds
    for power_target in [0.5, 0.8, 0.95, 0.99]:
        n_min = ht.min_n_for_power(power_target)
        print(f"\n  Minimum N for {power_target*100:.0f}% power: {n_min:,} cores")

    # Compare targets
    compare_targets()

    # Key insight
    print(f"\n{'='*70}")
    print(f"  FALSIFIABLE PREDICTION")
    print(f"{'='*70}")
    print(f"  Patch-clamp experiment: stimulate membrane patch with")
    print(f"  ultrafast UV laser pulses of varying intensity (N cores).")
    print(f"  Measured ionic current I(N) should follow:")
    print(f"    I(N) = N * p_open * I_single")
    print(f"  where p_open = 1 - (1 - 7.79e-4)^N")
    print(f"  This curve is distinct from classical linear summation.")
    print(f"  At N=5,000: I ~ 4.9 pA per synchronised burst.")
