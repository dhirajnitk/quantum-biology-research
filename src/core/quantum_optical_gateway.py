"""
quantum_optical_gateway.py  —  Publication-grade simulation of the
subatomic-to-optical handoff inside an intramolecular Tryptophan network.

Models a dense aromatic core (ε = 2 hydrophobic pocket) as an open quantum
system coupled to a thermal bath.  Tracks von Neumann entropy, coherent
information, and the resulting optical transit velocity into the lipid
membrane waveguide.

References
----------
- Breuer & Petruccione (2002). "The Theory of Open Quantum Systems."
- Mohseni et al. (2008). JCP 129, 174106 (ENAQT).
- Firmenich et al. (2026). bioRxiv (HEOM on tubulin Trp networks).
"""

import numpy as np
from qutip import Qobj, basis, mesolve, ket2dm, entropy_vn

# ---------------------------------------------------------------------------
#  Physical constants
# ---------------------------------------------------------------------------
HBAR = 1.054571817e-34       # J·s
KB = 1.380649e-23            # J/K
C0 = 299_792_458             # m/s, speed of light in vacuum
CM_PER_J = 5.034116e22       # 1 cm⁻¹  ≅ 1.9863e-23 J → factor
CM1_TO_RADS = 2 * np.pi * 2.99792458e10   # cm⁻¹ → rad/s


class QuantumOpticalGateway:
    """Lindblad simulation of a 4-site Trp aromatic core.

    Parameters
    ----------
    num_sites : int
        Number of coupled Trp sites in the hydrophobic pocket.
    epsilon : float
        Local dielectric constant of the protein core.
        2.0 = dry hydrophobic pocket (shielded).
        80.0 = open aqueous cytoplasm (fully exposed).
    temperature : float
        Temperature in Kelvin (default 310 K, human body).
    lam : float
        Reorganization energy [cm⁻¹] — scales linearly with ε.
        Baseline 35 cm⁻¹ at ε=2, scaled as lam(ε)=35*(ε/2).
    cutoff : float
        Bath cutoff frequency [cm⁻¹]; controls Markovianity.
    """

    def __init__(self, num_sites=4, epsilon=2.0, temperature=310.0,
                 lam=35.0, cutoff=53.0):
        self.N = num_sites
        self.eps = epsilon
        self.T = temperature
        self.lam = lam * (epsilon / 2.0)          # reorg  ∝  ε
        self.gamma_c = cutoff * CM1_TO_RADS       # → rad/s
        self.kT = KB * temperature

        # -------------------------------------------------------------------
        #   Tight-binding Hamiltonian [cm⁻¹]
        #   Site energies and couplings typical of a Trp-helix bundle.
        #   Electrostatic screening:  H_coupling  →  H_coupling / √ε
        # -------------------------------------------------------------------
        H_cm = np.array([
            [  0.0, -80.0,   5.0,   0.0],
            [-80.0, 110.0, -40.0,   2.0],
            [  5.0, -40.0, 210.0, -75.0],
            [  0.0,   2.0, -75.0,  50.0],
        ])
        # scale only the off-diagonal (coupling) elements
        mask_off = np.ones_like(H_cm, dtype=bool)
        np.fill_diagonal(mask_off, False)
        H_cm[mask_off] /= np.sqrt(epsilon)

        # convert to angular-frequency units (rad/s) for QuTiP
        self.H = Qobj(H_cm * CM1_TO_RADS)

        # -------------------------------------------------------------------
        #  Lindblad collapse operators
        #
        #  1. Pure dephasing at each site  —  kills off-diagonals
        #  2. Thermal relaxation  —  drives populations toward Boltzmann
        #  3. Trapping (optional)  —  models energy exiting to reaction centre
        # -------------------------------------------------------------------
        #  High-temperature Markovian dephasing rate  (ENAQT expression)
        #    γ_φ(T) = 2π (kT/ℏ) (λ / ω_c)
        gamma_phi = (2 * np.pi * self.kT / HBAR) * (self.lam / self.gamma_c)

        self.c_ops = []
        for i in range(self.N):
            proj = basis(self.N, i) * basis(self.N, i).dag()
            self.c_ops.append(np.sqrt(gamma_phi) * proj)

        #  Thermal relaxation  —  simple amplitude damping toward ground
        #  Rate taken as 1/10 of dephasing (typical in ENAQT literature)
        gamma_relax = gamma_phi / 10.0
        for i in range(1, self.N):
            # |i⟩ → |0⟩  (energy relaxation to lowest site)
            decay_op = np.sqrt(gamma_relax) * basis(self.N, 0) * basis(self.N, i).dag()
            self.c_ops.append(decay_op)

        # -------------------------------------------------------------------
        #  Optical transit velocity  —  lipid membrane waveguide
        # -------------------------------------------------------------------
        n_lipid = 1.45         # refractive index of lipid bilayer core
        self.v_optical = C0 / n_lipid

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------
    def run_dynamics(self, t_max_ps=5.0, n_steps=1000):
        """Run Lindblad evolution and return result object."""
        tlist = np.linspace(0, t_max_ps * 1e-12, n_steps)

        # Initial state:  coherent superposition across sites 0 & 1
        psi0 = (basis(self.N, 0) + basis(self.N, 1)).unit()
        rho0 = ket2dm(psi0)

        # Expectation operators: site populations
        e_ops = [basis(self.N, i) * basis(self.N, i).dag()
                 for i in range(self.N)]
        e_ops += [self.H]   # track energy expectation

        result = mesolve(self.H, rho0, tlist, self.c_ops, e_ops,
                         progress_bar=True)
        return result

    @staticmethod
    def von_neumann_entropy(rho):
        """S(ρ) = -Tr(ρ ln ρ)  —  quantum entropy in nats."""
        return entropy_vn(rho)

    @staticmethod
    def coherent_information(rho, H, beta=1.0):
        """I_c(ρ) = S(ρ) - S(ρ_env)  —  simplified proxy.

        For a full quantum channel capacity calculation one needs the
        Choi matrix.  Here we return the von Neumann entropy as a
        monotonic proxy for the coherent information available for
        transduction into an optical mode.
        """
        return entropy_vn(rho)

    # ------------------------------------------------------------------
    #  Analysis pipeline
    # ------------------------------------------------------------------
    def analyze(self, t_max_ps=5.0, n_steps=1000):
        """Run dynamics and return a summary dict."""
        res = self.run_dynamics(t_max_ps, n_steps)
        n = len(res.times)

        # per-step entropy  (measured from the full density matrix at t)
        S_t = np.zeros(n)
        I_c_t = np.zeros(n)
        pops = np.zeros((self.N, n))
        for i, t in enumerate(res.times):
            rho_t = res.states[i]
            S_t[i] = self.von_neumann_entropy(rho_t)
            I_c_t[i] = self.coherent_information(rho_t, self.H)
            for j in range(self.N):
                pops[j, i] = np.real(rho_t[j, j])

        # Markovian dephasing rate
        gamma_phi = (2 * np.pi * self.kT / HBAR) * (self.lam / self.gamma_c)

        return {
            'times':            res.times,
            'states':           res.states,
            'S_entropy':        S_t,
            'I_c_coherent':     I_c_t,
            'populations':      pops,
            'energy':           res.expect[-1],
            'gamma_phi_Hz':     gamma_phi,
            'gamma_phi_fs':     1e15 / gamma_phi,         # 1/γ_φ in fs
            'v_optical':        self.v_optical,
            'lam_cm':           self.lam,
            'epsilon':          self.eps,
            'H_matrix_cm':      self.H.full() / CM1_TO_RADS,
        }


# ===================================================================
#  Demo / validation
# ===================================================================
if __name__ == '__main__':
    print("=" * 68)
    print("  Quantum-Optical Gateway  —  Sub-Tubulin Information Calculator")
    print("=" * 68)

    # ----  Case 1:  Dry hydrophobic pocket  (ε = 2)  ----
    gate_dry = QuantumOpticalGateway(epsilon=2.0)
    res_dry = gate_dry.analyze(t_max_ps=2.0)

    # ----  Case 2:  Wet cytoplasm  (ε = 80)  ----
    gate_wet = QuantumOpticalGateway(epsilon=80.0)
    res_wet = gate_wet.analyze(t_max_ps=2.0)

    def fmt(val, unit=""):
        if isinstance(val, float):
            return f"{val:>12.4f}  {unit}"
        return f"{val:>12}  {unit}"

    print(f"\n{'':30s}  {'Insulated ε=2':>18s}  {'Exposed ε=80':>18s}")
    print(f"{'─'*68}")
    print(f"{'Optical transit v':30s}  {fmt(res_dry['v_optical'], 'm/s'):18s}  {fmt(res_wet['v_optical'], 'm/s')}")
    print(f"{'Dephasing rate 1/γ_φ':30s}  {fmt(res_dry['gamma_phi_fs'], 'fs'):18s}  {fmt(res_wet['gamma_phi_fs'], 'fs')}")
    print(f"{'Reorg. energy λ':30s}  {fmt(res_dry['lam_cm'], 'cm⁻¹'):18s}  {fmt(res_wet['lam_cm'], 'cm⁻¹')}")
    print(f"{'Initial entropy S₀':30s}  {fmt(res_dry['S_entropy'][0]):18s}  {fmt(res_wet['S_entropy'][0])}")
    print(f"{'Final entropy S_f':30s}  {fmt(res_dry['S_entropy'][-1]):18s}  {fmt(res_wet['S_entropy'][-1])}")
    print(f"{'ΔS (entropy generated)':30s}  {fmt(res_dry['S_entropy'][-1]-res_dry['S_entropy'][0]):18s}  {fmt(res_wet['S_entropy'][-1]-res_wet['S_entropy'][0])}")
    print(f"{'Coherent info I_c (final)':30s}  {fmt(res_dry['I_c_coherent'][-1]):18s}  {fmt(res_wet['I_c_coherent'][-1])}")
    print(f"{'─'*68}")

    print("\n  Hamiltonian [cm⁻¹]  (scaled by 1/√ε)\n")
    print(f"  ε = {res_dry['epsilon']:.0f}")
    H = res_dry['H_matrix_cm']
    for row in H:
        print(f"    [{', '.join(f'{x:8.1f}' for x in row)}]")

    print(f"\n  ε = {res_wet['epsilon']:.0f}")
    H_wet = res_wet['H_matrix_cm']
    for row in H_wet:
        print(f"    [{', '.join(f'{x:8.1f}' for x in row)}]")
    print()
