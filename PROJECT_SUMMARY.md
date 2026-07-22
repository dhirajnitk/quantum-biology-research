# Project Summary: Sub-Tubular Quantum Information Processing

## A Multi-Scale Architecture for Neural Computation via Tryptophan Networks

**Author:** Dhiraj Kumar
**Repository:** https://github.com/dhirajnitk/subnanometer-trp-transport
**Latest commit:** `e60bc7d` (27 files, 2985 insertions)
**Date:** July 2026

---

## Table of Contents

1. [The Hypothesis](#1-the-hypothesis)
2. [Repository Structure](#2-repository-structure)
3. [Code Modules — Complete Reference](#3-code-modules--complete-reference)
4. [Key Findings — All Numerical Results](#4-key-findings--all-numerical-results)
5. [Bug Fixes During Development](#5-bug-fixes-during-development)
6. [Paper Roadmap](#6-paper-roadmap)
7. [How to Run Everything](#7-how-to-run-everything)
8. [References](#8-references)

---

## 1. The Hypothesis

### 1.1 The 20-Watt Paradox

The human brain consumes approximately **20 W** of metabolic power, yet performs massively parallel computations — high-fidelity visual interpretation, motor coordination, memory retrieval, abstract thought — that would require **megawatts** in a silicon supercomputer. This efficiency gap of **six orders of magnitude** cannot be explained by classical action potentials alone, which are slow (1-100 m/s) and energetically expensive (~2.4 × 10⁵ ATP molecules per action potential).

### 1.2 Why Microtubule Models Fail

Previous proposals for quantum effects in the brain focused on microtubule networks (Orch-OR theory). However:
- **Tegmark (1996):** Calculated decoherence time of ~10⁻¹³ s for microtubule superpositions
- **Firmenich et al. (2026):** Rigorous HEOM calculations show ~13 fs dephasing at 310 K for the 1JFF tubulin Trp network
- **Energy budget:** GTP hydrolysis falls ~10⁵ short of the power needed for macro-scale coherence

### 1.3 The Sub-Tubular Solution

Quantum information processing occurs not in microtubule lattices, but in **sub-tubular hydrophobic pockets** within individual synaptic proteins. Dense Tryptophan (Trp) aromatic networks, shielded by the low-dielectric protein core (ε ≈ 2), maintain coherence for 50-100 fs — sufficient for ultra-fast local computation but not for long-range quantum transmission.

**The tri-layer architecture:**

```
Layer 1 (Quantum):   Trp networks in ε≈2 hydrophobic pockets
    |                50-100 fs coherence, exciton hopping < 1.5 nm
    |                Biophoton emission at 327 nm on collapse
    ▼
Layer 2 (Optical):   Lipid membrane waveguide (n=1.45)
    |                Total internal reflection, speed-of-light relay
    |                Z-channel with zero false positives
    ▼
Layer 3 (Electrical): Ion channel gating via allosteric trigger
                     Macro-scale action potentials
```

**The paradigm shift:** Free-space energy transfer between Trp cores is physically impossible (p ≈ 1.3 × 10⁻⁵ per exciton). The biophoton acts as a **binary information trigger** in a Z-channel, not as an energy source. Spatial ensemble of 5,000-10,000 parallel cores achieves >94% gating reliability.

---

## 2. Repository Structure

```
quantum-biology-research/
│
├── docs/
│   ├── reference/
│   │   └── QUANTUM_BIOLOGY_REFERENCE.md   # 2500+ line reference (Parts I-XI)
│   ├── architecture/
│   │   └── Architecture_of_Information_...md  # Narrative essay
│   └── research_plan/
│       └── RESEARCH_PLAN.md              # 12-section research plan
│
├── src/
│   ├── pdb_tools/
│   │   ├── trp_extractor.py             # PDB download + Trp extraction + distances
│   │   └── batch_processor.py           # Batch analysis + KCKAS + A_crit
│   │
│   ├── core/
│   │   ├── hamiltonian.py               # H_total = H_PDB + H_clock(A) builder
│   │   ├── lindblad_solver.py           # FMO 7-site Liouville Lindblad solver
│   │   ├── biophoton_relay.py           # Z-channel + spatial ensemble relay
│   │   ├── quantum_optical_gateway.py   # QuTiP Lindblad (QuTiP DLL blocked)
│   │   ├── quantum_optical_gateway_simple.py  # NumPy pedagogical version
│   │   └── neuro_registry.py            # PDB target metadata dictionary
│   │
│   ├── analysis/
│   │   ├── pdb_contextuality.py         # KCKAS contextuality on real PDB geometry
│   │   ├── quantum_hamiltonian_engine.py # Clock-driven coherence rescue + S(A)
│   │   ├── phase_pumped_kckas.py        # DEPRECATED — use HamiltonianEngine instead
│   │   ├── pdb_bell_test.py             # CHSH Bell inequality (encoding artifacts)
│   │   ├── z_channel_capacity.py        # Shannon capacity, Landauer energy
│   │   ├── hypothesis_test.py           # Statistical gating framework
│   │   ├── quantum_signature_test.py    # Superradiance/phase/CHSH tests
│   │   ├── manuscript_compiler.py       # Unified data matrix for publication
│   │   ├── channel_capacity.py          # Holevo + entropy production (P3)
│   │   ├── quantum_darwinism.py         # Pointer states + redundancy + SBS (P2)
│   │   ├── p3_thermodynamics.py         # Full P3 standalone analysis
│   │   └── p1_fmo/
│   │       └── fmo_lindblad.py          # FMO 7-site ENAQT + QMI
│   │
│   ├── ml/
│   │   ├── gnn_transport_predictor.py   # StructuralQuantumGnn + training (P4)
│   │   └── gnn_pipeline.py             # Simplified GCN for single target
│   │
│   └── main_fmo_pipeline.py            # Master integration script (P1-P4)
│
├── papers/
│   ├── p0_sub_tubulin/
│   │   └── outline.md                  # P0 paper outline
│   └── sub_tubulin_manuscript/
│       ├── manuscript.tex              # Physical Review E manuscript draft
│       ├── generate_figures.py         # Generates 4 manuscript figures
│       └── figures/
│           ├── fig1_z_channel.pdf       # Z-channel diagram + ensemble curve
│           ├── fig2_geometry_tax.pdf    # KCKAS S vs symmetry bar chart
│           ├── fig3_hamiltonian.pdf     # S(A) curves for 6 targets
│           └── fig4_energy.pdf          # Landauer efficiency comparison
│
├── data/
│   ├── cache/                          # PDB structure files (gitignored)
│   ├── 7TYO_trp_network.csv           # NMDA receptor Trp distances
│   ├── 1BL8_trp_network.csv           # KcsA K+ channel Trp distances
│   └── batch_results.csv              # Batch analysis output
│
├── PROJECT_SUMMARY.md                 # This file
├── RESEARCH_SUMMARY.md                 # 6-section concise summary
├── PROJECT_PLAN.md                    # Development roadmap
├── RESEARCH_PLAN.md                   # 12-section research plan
├── README.md                          # Repository overview
├── requirements.txt                   # Dependencies (numpy, scipy, qutip, etc.)
└── LICENSE                            # MIT
```

---

## 3. Code Modules — Complete Reference

### 3.1 PDB Tools

#### `src/pdb_tools/trp_extractor.py`
- Downloads PDB files from RCSB (with local cache)
- Extracts Trp CG atom coordinates
- Computes pair-distance matrix
- Classifies pairs: "coupled" (< 1.5 nm) or "optical relay" (1.5-5 nm)
- Builds tight-binding Hamiltonian from distances: `J_ij = J0 × (R0/R_ij)³`
- **Usage:** `python src/pdb_tools/trp_extractor.py 7TYO --chain A --save`

#### `src/pdb_tools/batch_processor.py`
- Production batch analysis using HamiltonianEngine
- For each target: download → extract → compute distance matrix → build Hamiltonian → KCKAS → find A_crit
- Outputs formatted manuscript table with all metrics
- **13 targets analysed:** 6PV7, 1BL8, 7TYO, 6LQA, 6CNO, 7KOX, 6J8J, 1YAG, 1JFF, 3N2K, 5VKH, 4UVN, 3JBR

### 3.2 Core Modules

#### `src/core/hamiltonian.py`
- Formal multi-scale Hamiltonian builder
- `build_pdb_hamiltonian(xyz_coords, dielectric)` — Dexter-type couplings scaled by 1/√ε
- `build_clock_hamiltonian(distances, amplitude)` — exponential phase pump: J = A × 400 × exp(-R/15)
- `build_multiscale_hamiltonian(...)` — H_total = H_PDB + H_clock(A)
- `coherence_factor(H_cm)` — F_coh = Σ|H_ij| / (Σ|H_ij| + Σ|H_ii|)
- `kckas_s(H_cm)` — S = F_coh × √5

#### `src/core/lindblad_solver.py`
- FMO 7-site Lindblad solver with Liouville-space propagation
- Adolphs & Renger (2006) Hamiltonian
- Non-Hermitian trapping at site 3: H_eff = H - i × k_trap/2 × |2⟩⟨2|
- Pure dephasing via |k⟩⟨k| Lindblad operators
- Static disorder parameter (cm⁻¹) for Anderson localization
- `run_time_evolution(dephasing_rate, t_max, dt)` → efficiency, density_history
- `compute_qmi_matrix(rho)` → 7×7 QMI matrix
- `enaqt_curve(gamma_range)` → ENAQT efficiency vs dephasing

**Critical bug fix:** Liouvillian must use `H.conj()` (complex conjugate) not `H.T` (transpose) for non-Hermitian H. Using `H.T` gave incorrect dynamics where trap population never decayed (trace always stayed at 1.0).

#### `src/core/biophoton_relay.py`
- Z-channel model with spatial ensemble summation
- Corrected physics: QY = 0.21 (measured, not derived), exciton = 35,700 cm⁻¹ (280 nm)
- Trp target: p = 1.34 × 10⁻⁵ per core
- CCO target: p = 5.92 × 10⁻⁴ per core
- `spatial_ensemble_success(N_cores, target_type)` → P_ensemble, p_hit, n_arr

### 3.3 Analysis Modules

#### `src/analysis/pdb_contextuality.py`
- KCKAS contextuality test using 5 Trp residues forming a pentagram
- Builds optimal 5-vector projectors on 3D Hilbert space
- Coherence factor from coupling-to-disorder ratio
- S = F_coh × √5. Classical bound: S ≤ 2. Quantum max: S ≤ 2.236
- **Note:** This is a "KCKAS coherence score" — a heuristic based on geometric coupling, not a full quantum state optimisation

#### `src/analysis/quantum_hamiltonian_engine.py`
- Clock-driven coherence rescue
- `compute_kckas_s(A)` — exponential blend: S(A) = S₀ + (√5 - S₀) × (1 - exp(-J_eff/ΔE))
- `scan_clock_drive(drives)` → S values at each amplitude
- `find_critical_drive()` → binary search for minimum A achieving S > 2
- **Key result:** All 8 PDB targets breach S > 2 under full clock drive

#### `src/analysis/channel_capacity.py`
- `FmoThermodynamicEngine` class
- `compute_holevo_capacity(time_index)` — HSW capacity χ
- `compute_entropy_production(time_index)` — dS/dt = -Tr(dρ/dt × log ρ)
- `compute_entropy_trajectory()` — S(t)
- `compute_energy_per_bit(time_index, efficiency)` — Landauer cost
- `classical_random_walk_efficiency()` — FRET-like comparison

#### `src/analysis/quantum_darwinism.py`
- `QuantumDarwinism` class (also aliased as `FmoDarwinismAnalyzer`)
- `compute_redundancy_curve(time_index)` — partial information curve
- `redundancy_at_delta(delta, time_index)` — R(δ)
- `find_pointer_states(time_index)` — einselection via env eigendecomposition
- `compute_sbs_fidelity(time_index)` — spectrum broadcast structure
- `mi_between_sites(i, j, time_index)` — pairwise quantum mutual information

### 3.4 Machine Learning

#### `src/ml/gnn_transport_predictor.py`
- `StructuralQuantumGnn` — 2-layer GCN with backpropagation
- Input features: [site energy, dielectric, SASA proxy, connectivity degree]
- Dual output: [ENAQT efficiency, Holevo capacity]
- `train_step(A, X, targets, lr)` — full analytical gradient backprop
- `QuantumTrainingDataFactory` — generates synthetic training data
- **Result:** Loss drops from 178 → 0.036 over 500 epochs; test error ~10%

### 3.5 Integration

#### `src/main_fmo_pipeline.py`
- Runs all 4 pipeline stages sequentially:
  1. Lindblad solver (ENAQT + QMI)
  2. Quantum Darwinism (PIC + pointer states + SBS)
  3. Thermodynamics (Holevo + entropy production)
  4. GNN prediction
- Outputs manuscript-ready data rows

---

## 4. Key Findings — All Numerical Results

### 4.1 Structural Trp Networks in Neural Receptors

| PDB | Protein | Class | Trp | Coupled | Relay | Dielectric | Trp pair data file |
|-----|---------|-------|:---:|:-------:|:-----:|:----------:|:------------------:|
| 1BL8 | KcsA K+ channel | V-gated | 5 | 3 | 7 | 2.0 | `data/1BL8_trp_network.csv` |
| 6PV7 | nACh receptor | L-gated | 26 | 22 | 303 | 2.0 | |
| 7TYO | NMDA receptor GluN1/GluN2A | L-gated | 23 | 21 | 232 | 2.1 | `data/7TYO_trp_network.csv` |
| 6LQA | NaV1.4 sodium channel | V-gated | 25 | 18 | 282 | 2.3 | |
| 6CNO | NMDA receptor GluN1/GluN2B | L-gated | 10 | 5 | 40 | 2.1 | |
| 7KOX | Alpha-7 nicotinic AChR | L-gated | 12 | 4 | 62 | 2.0 | |
| 6J8J | NaV1.4 (alternative) | V-gated | 27 | 20 | 331 | 2.3 | |
| 1YAG | SNARE complex | Synaptic | 6 | 3 | 12 | 2.6 | |
| 1JFF | Tubulin dimer (bovine) | Cytoskeletal | 5 | 1 | 9 | 2.1 | |
| 3N2K | Tubulin dimer (mammalian) | Cytoskeletal | 5 | 1 | 9 | 2.1 | |

### 4.2 Z-Channel Gating — Spatial Ensemble Summation

#### Corrected physics parameters
| Parameter | Old (wrong) | New (correct) | Source |
|-----------|:-----------:|:-------------:|--------|
| Trp quantum yield | 0.82 (sqrt(80/ε) × 0.13) | **0.21** | Babcock (2024), measured |
| Exciton energy | 12,000 cm⁻¹ (830 nm) | **35,700 cm⁻¹** (280 nm) | Babcock (2024), Trp absorption peak |
| Trp target p_per_core | 1.80 × 10⁻⁵ | **1.34 × 10⁻⁵** | Corrected QY + exciton |
| CCO target p_per_core | 7.79 × 10⁻⁴ | **5.92 × 10⁻⁴** | Corrected QY + exciton |

#### Ensemble gating with CCO target

| N_cores | P_success | Channel capacity | N for threshold |
|:-------:|:---------:|:---------------:|:--------------:|
| 1 | 0.06% | 0.993 bits | — |
| 100 | 5.7% | 0.764 bits | — |
| 1,000 | 44.7% | 0.010 bits | N₅₀ ≈ 1,173 |
| **5,000** | **94.8%** | **0.707 bits** | N₉₅ ≈ 5,060 |
| **10,000** | **99.7%** | **0.971 bits** | N₉₉ ≈ 7,780 |
| 25,000 | 99.999% | 0.999 bits | — |

#### Energy efficiency
| Metric | Value |
|--------|:-----:|
| Energy per bit at N=5,000 | 1.26 × 10⁻¹⁵ J |
| Landauer ratio | **4.1 × 10⁵ ×** |
| Classical CMOS | ~10⁷-10⁸ × Landauer |
| **Z-channel advantage** | **~100× more efficient than CMOS** |

### 4.3 KCKAS Contextuality — Static Geometry

| Target | Static S | F_coh | Mean distance (Å) | Coupled pairs |
|--------|:-------:|:-----:|:-----------------:|:-------------:|
| 1BL8 (KcsA) | **1.9735** | 0.99 | 25.4 | 3 |
| 6PV7 (nAChR) | 1.9174 | 0.96 | 59.0 | 22 |
| 7TYO (NMDA) | 1.8092 | 0.90 | 53.5 | 21 |
| 6J8J (NaV) | 1.7323 | 0.87 | 52.3 | 20 |
| 6LQA (NaV) | 1.7236 | 0.86 | 47.6 | 18 |
| 7KOX (nAChR) | 1.5824 | 0.79 | 48.1 | 4 |
| 6CNO (NMDA) | 1.3504 | 0.68 | 41.1 | 5 |
| 1YAG (SNARE) | 1.3124 | 0.66 | 25.5 | 3 |

**Classical bound: S ≤ 2.0 — no static structure breaches it.**

### 4.4 Hamiltonian Engine — Clock-Driven Coherence Rescue

| Target | S_static | A=0.25 | A=0.50 | A=0.75 | **A=1.00** | **A_crit** |
|--------|:-------:|:------:|:------:|:------:|:---------:|:----------:|
| 1BL8 (KcsA) | 1.97 | **2.12** | **2.19** | **2.22** | **2.23** | **0.03** |
| 6PV7 (nAChR) | 1.92 | **2.02** | **2.09** | **2.13** | **2.17** | **0.20** |
| 7TYO (NMDA) | 1.81 | 1.95 | **2.05** | **2.11** | **2.15** | **0.36** |
| 6LQA (NaV) | 1.72 | 1.87 | 1.98 | **2.05** | **2.11** | **0.57** |
| 6CNO (NMDA) | 1.35 | 1.68 | 1.89 | **2.02** | **2.10** | **0.71** |
| 7KOX (nAChR) | 1.58 | 1.81 | 1.96 | **2.05** | **2.12** | **0.61** |
| 6J8J (NaV) | 1.73 | 1.85 | 1.94 | **2.01** | **2.06** | **0.71** |
| 1YAG (SNARE) | 1.31 | 1.70 | 1.92 | **2.05** | **2.13** | **0.63** |

**All targets breach S > 2.0 at full clock drive (A = 1).**

### 4.5 FMO ENAQT — Thermodynamics

#### ENAQT with static disorder (50 cm⁻¹)

| γ (cm⁻¹) | γ (ps⁻¹) | Efficiency | Regime |
|:--------:|:--------:|:----------:|--------|
| 0.1 | 0.002 | 9.0% | Anderson localization |
| 10 | 0.19 | 15.4% | Coherent |
| 50 | 0.94 | 26.9% | ENAQT ascent |
| **100** | **1.88** | **29.4%** | **ENAQT peak (3.27× enhancement)** |
| 175 | 3.29 | 27.5% | ENAQT descent |
| 500 | 9.40 | 17.0% | Quantum Zeno |

#### P3 Thermodynamic comparison

| γ (cm⁻¹) | Efficiency | Peak S | Peak χ | Peak dS/dt | × Landauer |
|:--------:|:----------:|:------:|:------:|:----------:|:----------:|
| 0 | 8.9% | 0.12 | 0.12 | 0.000 | 2.26 × 10⁴ |
| 10 | 15.4% | 2.07 | 2.07 | 0.247 | 775× |
| 50 | 26.9% | 2.26 | 2.26 | 1.050 | 406× |
| **100** | **29.4%** | 2.17 | 2.17 | 1.470 | **387×** |
| 175 | 27.5% | 2.07 | 2.07 | 1.709 | 434× |
| 500 | 17.0% | 1.81 | 1.81 | 1.639 | 799× |

### 4.6 Quantum Darwinism (FMO at γ = 175 cm⁻¹)

#### Partial Information Curves

| Fragment size | I(sys:frag) at t=0.5ps | I(sys:frag) at t=2.0ps | I(sys:frag) at t=10ps |
|:-------------:|:----------------------:|:----------------------:|:---------------------:|
| 1 | 0.0525 | 0.0078 | 0.0008 |
| 2 | 0.0705 | 0.0148 | 0.0017 |
| 3 | 0.0885 | 0.0192 | 0.0022 |
| 4 | 0.1155 | 0.0256 | 0.0030 |
| 5 | 0.1594 | 0.0355 | 0.0043 |
| 6 | 0.2716 | 0.0427 | 0.0047 |

#### Pairwise QMI (site 1 with each environment site)

| Pair | I (bits) |
|:----:|:--------:|
| Site 1 - Site 2 | 0.7311 |
| Site 1 - Site 6 | 0.4821 |
| Site 1 - Site 2 | 0.7311 |
| Site 1 - Site 7 | 0.4487 |
| Site 1 - Site 3 | 0.4201 |
| Site 1 - Site 5 | 0.4298 |
| Site 1 - Site 4 | 0.4341 |

**SBS fidelity: 1.0 — environment fully encodes system information.**

#### Pointer states
| λ₁ | λ₂ | λ₃ | λ₄ | λ₅ |
|:--:|:--:|:--:|:--:|:--:|
| 0.310 | 0.080 | 0.053 | 0.043 | 0.040 |

### 4.7 GNN Transport Prediction (P4)

| Metric | Value |
|--------|:-----:|
| Training epochs | 500 |
| Initial loss | 178.5 |
| Final loss | 0.036 |
| Test ENAQT error | 10.2% |
| Test Holevo error | 13.4% |
| Architecture | 2-layer GCN (4→12→2) |
| Training data | Synthetic (from FMO parameter space) |

### 4.8 Geometry Tax — Summary

```
A_crit ranking:
  0.03  (1BL8, KcsA)    — Perfect homotetramer, tight ring geometry
  0.09  (5VKH)            — Similarly symmetric structure
  0.20  (6PV7, nAChR)    — Pentameric, bulky extracellular domain
  0.36  (7TYO, NMDA)     — Multi-domain complex, 40-60 Å gaps
  0.57  (6LQA, NaV)      — One polypeptide chain, 4 non-identical domains
  0.61  (7KOX, nAChR)    — Alpha-7 receptor
  0.63  (1YAG, SNARE)    — Linear bundle
  0.71  (6CNO, NMDA)     — Asymmetric hinge regions
  0.71  (6J8J, NaV)      — Alternative NaV conformation
  0.85  (1JFF, tubulin)  — Sparse Trp network, wide spacing
```

---

## 5. Paper Roadmap

| Paper | Title | Target Journal | Core Contribution | Status |
|-------|-------|---------------|-------------------|--------|
| **P0** | Sub-Tubular Quantum Information Processing: A Multi-Scale Architecture for Neural Computation via Tryptophan Networks | *Physical Review E* | Tri-layer architecture, Z-channel, KCKAS on 10 PDB targets, Hamiltonian engine, geometry tax | ✅ Manuscript at `papers/sub_tubulin_manuscript/manuscript.tex` (4 figures, 4 tables, 24 refs) |
| **P1** | Quantum Mutual Information Reveals Energy Transfer Pathways in the FMO Complex | *Journal of Chemical Physics* | ENAQT efficiency curve, QMI matrix across 21 chromophore pairs, optimal dephasing | ✅ Code at `src/core/lindblad_solver.py` + `src/analysis/p1_fmo/fmo_lindblad.py` |
| **P2** | Quantum Darwinism in Photosynthetic Energy Transfer | *New Journal of Physics* | Pointer states, redundancy R(δ), SBS fidelity, partial information curves | ✅ Code at `src/analysis/quantum_darwinism.py` |
| **P3** | Thermodynamic Cost of Quantum Coherence in Photosynthesis | *Physical Review E* or *PRX Life* | Holevo capacity, entropy production, Landauer cost, quantum vs classical comparison | ✅ Code at `src/analysis/channel_capacity.py` + `src/analysis/p3_thermodynamics.py` |
| **P4** | Machine Learning Prediction of Quantum Transport from Protein Structure | *npj Quantum Information* | GCN surrogate model, ENAQT + Holevo prediction from structural features | ✅ Code at `src/ml/gnn_transport_predictor.py` |

---

## 6. How to Run Everything

### Setup
```bash
# The project uses the Python environment at:
# C:\dhiraj\matrimony-dating-app\matrimony-platform\.venv_win\Scripts\python.exe

# Install dependencies
pip install -r requirements.txt
```

### P0: Core Hypothesis
```bash
# Extract Trp networks from a PDB structure
python src/pdb_tools/trp_extractor.py 7TYO --save

# Run batch analysis across all targets
python src/pdb_tools/batch_processor.py

# Run Hamiltonian engine with clock drive
python src/analysis/quantum_hamiltonian_engine.py

# Generate manuscript figures
python papers/sub_tubulin_manuscript/generate_figures.py

# Compile manuscript (requires LaTeX)
cd papers/sub_tubulin_manuscript && pdflatex manuscript.tex && pdflatex manuscript.tex
```

### P1-P4: FMO Pipeline
```bash
# Run the full FMO pipeline (all 4 papers)
python src/main_fmo_pipeline.py

# Run P3 thermodynamic analysis independently
python src/analysis/p3_thermodynamics.py

# Run P4 GNN transport predictor
python src/ml/gnn_transport_predictor.py
```

### Verification
```bash
# Check all modules import correctly
python -c "import sys; sys.path.insert(0,'src')
from core.hamiltonian import build_multiscale_hamiltonian
from core.lindblad_solver import FmoLindbladSolver
from analysis.quantum_darwinism import QuantumDarwinism, FmoDarwinismAnalyzer
from analysis.channel_capacity import FmoThermodynamicEngine
from ml.gnn_transport_predictor import StructuralQuantumGnn
print('All imports OK')"
```

---

## 7. References

1. Babcock et al. (2024) *J. Phys. Chem. B* **128**, 1525 — Trp fluorescence QY, spectra
2. Babcock et al. (2024) *Front. Phys.* **12**, 1234 — Superradiance in Trp neuroproteins
3. Gassab, Pusuluk & Craddock (2026) *Entropy* **28**, 204 — QI flow in Trp networks
4. Firmenich et al. (2026) *bioRxiv* — HEOM on 1JFF tubulin Trp network
5. Klyachko et al. (2008) *Phys. Rev. Lett.* **101**, 020403 — KCKAS contextuality
6. Clauser et al. (1969) *Phys. Rev. Lett.* **23**, 880 — CHSH Bell inequality
7. Adolphs & Renger (2006) *Biophys. J.* **91**, 2778 — FMO 7-site Hamiltonian
8. Mohseni et al. (2008) *J. Chem. Phys.* **129**, 174106 — ENAQT theory
9. Holevo (1998) *IEEE Trans. Info. Theory* **44**, 269 — HSW capacity
10. Landauer (1961) *IBM J. Res. Dev.* **5**, 183 — Landauer limit
11. Zurek (2003) *Rev. Mod. Phys.* **75**, 715 — Quantum Darwinism
12. van Meer et al. (2008) *Nat. Rev. Mol. Cell Biol.* **9**, 112 — Membrane optics
13. Patwa & Kurian (2026) *Phys. Rev. A* **113** — Superradiance in helical emitters
14. Fisher (2015) *Ann. Phys.* **362**, 593 — Posner molecules
15. Cover & Thomas (2006) *Elements of Information Theory*, 2nd ed. — Wiley
