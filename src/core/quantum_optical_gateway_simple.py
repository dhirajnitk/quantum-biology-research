"""
quantum_optical_gateway_simple.py  —  Pedagogical NumPy-only version.

Models a 4-site Trp aromatic core with a simple exponential-decay
approximation for decoherence.  Kept for conceptual clarity alongside
the full QuTiP Lindblad implementation in quantum_optical_gateway.py.
"""

import numpy as np


class QuantumOpticalGateway:
    def __init__(self, num_sites=4, dielectric_constant=2.0):
        self.num_sites = num_sites
        self.epsilon = dielectric_constant
        self.hbar = 1.0545718e-34
        self.kb = 1.380649e-23

        self.H_base = np.array([
            [  0.0, -80.0,   5.0,   0.0],
            [-80.0, 110.0, -40.0,   2.0],
            [  5.0, -40.0, 210.0, -75.0],
            [  0.0,   2.0, -75.0,  50.0]
        ])
        self.H = self.H_base / np.sqrt(self.epsilon)

    def compute_von_neumann_entropy(self, rho):
        evals = np.linalg.eigvalsh(rho)
        evals = evals[evals > 1e-15]
        return -np.sum(evals * np.log(evals))

    def run_phase_evolution(self, timesteps=100, dt_fs=1.0):
        rho = np.zeros((self.num_sites, self.num_sites), dtype=complex)
        rho[0, 0] = 0.5
        rho[1, 1] = 0.5
        rho[0, 1] = 0.5
        rho[1, 0] = 0.5

        entropy_ledger = []
        c_speed = 3e8 / np.sqrt(1.45)

        for t in range(timesteps):
            S = self.compute_von_neumann_entropy(rho)
            entropy_ledger.append(S)
            decay_factor = 0.999 - (0.005 * (self.epsilon / 2.0))
            rho = rho * decay_factor
            for i in range(self.num_sites):
                rho[i, i] = rho[i, i] / np.trace(rho)

        return entropy_ledger, c_speed


if __name__ == '__main__':
    for eps, label in [(2.0, "Insulated ε=2"), (80.0, "Exposed ε=80")]:
        gate = QuantumOpticalGateway(dielectric_constant=eps)
        S, v = gate.run_phase_evolution()
        print(f"[{label}]")
        print(f"  Optical transit velocity: {v:,.0f} m/s")
        print(f"  Initial entropy: {S[0]:.4f}")
        print(f"  Final entropy:   {S[-1]:.4f}")
        print(f"  Net entropy:     {S[-1] - S[0]:.4f}\n")
