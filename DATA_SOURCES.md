# Data Sources and Verification

This document tracks all empirical data sources, physical constants, and mathematical equations that ground the Sub-Tubular Quantum Information Processing framework. All numbers in the manuscript can be traced back to the sources below.

---

## 1. Structural Data — Protein Data Bank

| PDB ID | Protein | Class | Resolution | Trp residues | Source |
|--------|---------|-------|:----------:|:------------:|--------|
| 1BL8 | KcsA K⁺ channel | V-gated ion channel | 3.5 Å | 5 | rcsb.org |
| 6PV7 | Nicotinic ACh receptor | L-gated ion channel | 3.8 Å | 26 | rcsb.org |
| 7TYO | NMDA receptor GluN1/GluN2A | L-gated ion channel | 3.6 Å | 23 | rcsb.org |
| 6LQA | NaV1.4 sodium channel | V-gated ion channel | 3.2 Å | 25 | rcsb.org |
| 6CNO | NMDA receptor GluN1/GluN2B | L-gated ion channel | 3.3 Å | 10 | rcsb.org |
| 7KOX | α7 nicotinic AChR | L-gated ion channel | 3.1 Å | 12 | rcsb.org |
| 6J8J | NaV1.4 (alternative) | V-gated ion channel | 3.5 Å | 27 | rcsb.org |
| 1YAG | SNARE complex | Synaptic vesicle anchor | 2.4 Å | 6 | rcsb.org |
| 1JFF | Tubulin dimer (bovine) | Cytoskeletal | 3.5 Å | 5 | rcsb.org |
| 3N2K | Tubulin dimer (mammalian) | Cytoskeletal | 3.6 Å | 5 | rcsb.org |
| 5VKH | — | Membrane protein | — | 15 | rcsb.org |
| 4UVN | — | Membrane protein | — | — | rcsb.org |
| 3JBR | — | Membrane protein | — | 20 | rcsb.org |

**Extraction method:** CG (carbon gamma) atom coordinates as centroid proxy. Transition dipole vectors computed from heavy-atom coordinates of the indole ring (CG → CH2 direction for the S₀→S₁ transition long-axis approximation).

---

## 2. Optical and Spectroscopic Constants

| Parameter | Value | Source |
|-----------|:-----:|--------|
| Trp QY (in tubulin) | 0.21 | Babcock et al. (2024) JPCB |
| Trp absorption peak | 280 nm | Babcock et al. (2024) |
| Trp emission peak | 327 nm | Babcock et al. (2024) |
| Trp extinction coefficient at 280 nm | 5600 M⁻¹cm⁻¹ | Lakowicz (2006) |
| Trp absorption cross-section at 280 nm | 2.1 × 10⁻²¹ m² | Derived from ε₂₈₀ |
| Trp absorption cross-section at 327 nm | 1.8 × 10⁻²² m² | Gaussian lineshape from ε₂₈₀ |
| CCO absorption cross-section at 327 nm | 8.0 × 10⁻²¹ m² | Wikström et al. (2012) BBA |
| Lipid refractive index | 1.45 | van Meer et al. (2008) |
| Cytoplasm refractive index | 1.33 | Handbook of Biological Optics |
| Local dielectric (hydrophobic core) | 2.0–2.6 | Huang (2019) |

---

## 3. Physical Constants

| Constant | Value | Notes |
|----------|:-----:|-------|
| h̄ (reduced Planck) | 1.054571817 × 10⁻³⁴ J·s | CODATA |
| k_B (Boltzmann) | 1.380649 × 10⁻²³ J/K | CODATA |
| T (body temperature) | 310 K | Physiological |
| kT ln 2 (Landauer limit at 310 K) | 2.97 × 10⁻²¹ J | Landauer (1961) |
| Exciton energy (280 nm) | 35,700 cm⁻¹ = 7.09 × 10⁻¹⁹ J | Trp absorption peak |
| Photon energy (327 nm) | 30,600 cm⁻¹ = 6.08 × 10⁻¹⁹ J | Trp emission peak |
| c (speed of light) | 2.998 × 10⁸ m/s | CODATA |
| 1 cm⁻¹ → ps⁻¹ (h̄ = 1) | 0.0188 ps⁻¹ | Conversion factor |

---

## 4. Key Equations

### Dipole-dipole coupling (with orientation factor)
$$J_{ij} = \frac{J_0}{\sqrt{\varepsilon}} \left(\frac{R_0}{R_{ij}}\right)^3 \cdot |\kappa|, \quad \kappa = \hat{\mu}_i \cdot \hat{\mu}_j - 3(\hat{\mu}_i \cdot \hat{R}_{ij})(\hat{\mu}_j \cdot \hat{R}_{ij})$$

### Z-channel matrix
$$\mathbf{P} = \begin{pmatrix} 1 & 0 \\ 1-p & p \end{pmatrix}$$

### Spatial ensemble summation
$$P_{\text{success}}(N) = 1 - (1-p)^N$$

### KCKAS coherence score
$$S = F_{\mathrm{coh}} \times \sqrt{5}, \quad F_{\mathrm{coh}} = \frac{\sum_{i\neq j} |H_{ij}|}{\sum_{i\neq j} |H_{ij}| + \sum_i |H_{ii}|}$$

### Multi-scale Hamiltonian
$$H_{\mathrm{total}} = H_{\mathrm{PDB}} + H_{\mathrm{clock}}(A), \quad H_{\mathrm{clock}}(A) = \sum_{i<j} A \cdot J_{\max} e^{-R_{ij}/R_0}$$

### Clock-driven KCKAS rescue
$$S(A) = S_0 + (\sqrt{5} - S_0) \cdot \left[1 - \exp\!\left(-\frac{J_{\mathrm{eff}}(A)}{\Delta E}\right)\right]$$

### Energy per bit
$$E_{\mathrm{bit}} = \frac{N \times E_{\mathrm{exciton}}}{P_{\mathrm{success}}(N)}$$

### Whole-brain power
$$P_{\mathrm{Z}} = E_{\mathrm{bit}} \times N_{\mathrm{syn}} \times f$$

---

## 5. Verified Numerical Results

| Result | Value | How Verified |
|--------|:-----:|-------------|
| Trp-to-Trp re-excitation probability | 1.8 × 10⁻⁵ | Cross-section calculation from published optical data |
| CCO target per-core probability | 7.84 × 10⁻⁴ | Cross-section from Wikström (2012) |
| KcsA static KCKAS S | 1.9735 | PDB coordinate calculation |
| 1BL8 A_crit | 0.05 | HamiltonianEngine sweep |
| 6LQA A_crit | 0.58 | HamiltonianEngine sweep |
| ENAQT peak (50 cm⁻¹ disorder) | 29.4% at ~100 cm⁻¹ | Lindblad solver |
| Landauer ratio (neural) | 6.23 × 10⁵ × | Energy per bit calculation |
| Landauer ratio (FMO) | ~390 × | Channel capacity calculation |
| Whole-brain power (1 Hz) | 0.19 W | Extrapolation over 10¹⁴ synapses |
| Whole-brain power (10 Hz) | 1.9 W | Extrapolation over 10¹⁴ synapses |
| Maximum sustainable rate (10% budget) | 10.8 Hz | Power budget calculation |
| Total brain bandwidth | 0.87 Pbits/s | Channel capacity × synapses × rate |

---

## 6. Literature Sources

1. Babcock et al. (2024) *J. Phys. Chem. B* **128**, 1525 — Trp fluorescence QY and spectra
2. Lakowicz (2006) *Principles of Fluorescence Spectroscopy*, 3rd ed. — Trp extinction coefficients
3. Wikström et al. (2012) *Biochim. Biophys. Acta* **1817**, 77 — CCO absorption spectra
4. Landauer (1961) *IBM J. Res. Dev.* **5**, 183 — Landauer limit
5. Frank (2002) *Comput. Sci. Eng.* **4**, 16 — CMOS energy comparison
6. Adolphs & Renger (2006) *Biophys. J.* **91**, 2778 — FMO Hamiltonian
7. Klyachko et al. (2008) *PRL* **101**, 020403 — KCKAS contextuality
8. Firmenich et al. (2026) *bioRxiv* — HEOM on 1JFF tubulin Trp network
9. Mohseni et al. (2008) *JCP* **129**, 174106 — ENAQT theory
10. van Meer et al. (2008) *Nat. Rev. Mol. Cell Biol.* **9**, 112 — Membrane optics
11. Huang (2019) *Introduction to Statistical Physics* — Dielectric constants
12. Wong-Riley et al. (2005) *J. Biol. Chem.* **280**, 4761 — CCO photoreception
13. Karu (1999) *J. Photochem. Photobiol. B* **49**, 1 — CCO light absorption
14. Tang & Dai (2014) *Sci. Rep.* **4**, 5900 — Neural biophoton propagation
15. Dicke (1954) *Phys. Rev.* **93**, 99 — Superradiance theory
