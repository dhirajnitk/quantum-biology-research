# Project Plan: Information-Theoretic Protein Analysis & Sub-Tubulin Quantum Biology

---

## 1. Document Architecture — Final Structure

**QUANTUM_BIOLOGY_REFERENCE.md**

| Part | Title | Status |
|------|-------|--------|
| I | Foundations of Quantum Mechanics for Biology | ✅ Complete |
| II | Quantum Biology — The Field | ✅ Complete |
| III | The FMO Complex — The Gold Standard | ✅ Complete |
| IV | Sub-Tubulin Systems & Neural Quantum Effects | ✅ Complete |
| V | Other Quantum Biology Systems | ✅ Complete |
| VI | Computational Methods — Full Toolkit | 🟡 Needs update |
| VII | Information Theory & Quantum Biology | 🟡 Needs update |
| VIII | Researcher Profile — Dr. Anita Goel | ✅ Complete |
| IX | Research Gaps (2024-2026) | ✅ Complete |
| X | Annotated Bibliography | ✅ Complete |
| XI | Ongoing Investigations | 🟡 In progress |

**Section 11 sub-sections:**

| Sub-section | Title | Status |
|-------------|-------|--------|
| 11.1 | Multi-Modal Communication Architecture | ✅ Complete |
| 11.2 | Phase Synchronization for Complex Routing | ✅ Complete |
| 11.3 | Open Questions | ✅ Complete |
| 11.4 | The Virtual Lab: Implementation Roadmap | ✅ Complete |
| 11.5 | Geometrically Sub-Threshold Quantum Vault | ✅ Complete |
| 11.6 | Paper Roadmap | ✅ Complete |

**Other files in repository:**

| File | Role | Status |
|------|------|--------|
| `RESEARCH_PLAN.md` | Standalone 12-section research plan | ✅ Complete |
| `Architecture_of_Information_Quantum_Coherence_Cosmic_Entropy.md` | Narrative essay | ✅ Complete |
| `quantum_optical_gateway.py` | QuTiP Lindblad simulation (publication-grade) | ✅ Written |
| `quantum_optical_gateway_simple.py` | NumPy pedagogical version | ✅ Written |
| `PROJECT_PLAN.md` | This file | ✅ Written |

---

## 2. Paper-Writing Roadmap

### Target Journals by Tier

| Tier | Journal | Scope | OA |
|------|---------|-------|----|
| 1 | *Physical Review X* (PRX) | High-impact physics | Hybrid |
| 1 | *PRX Life* | New venue, well-aligned | Hybrid |
| 2 | *Journal of Chemical Physics* | Solid, reputable | Hybrid |
| 2 | *New Journal of Physics* | Open access, good fit | ✅ OA |
| 3 | *npj Quantum Information* | Nature portfolio, high impact | ✅ OA |
| 3 | *BioSystems* | Interdisciplinary, lower bar | Hybrid |

### Paper Pipeline

| Paper | Title | Target Journal | Core Tool | Timeline |
|-------|-------|---------------|-----------|----------|
| **P0** | Sub-Tubular Quantum Information Processing: A Multi-Modal Architecture for Neural Computation | *BioSystems* or *Quantum Biology* | Theoretical framework | Month 3-4 |
| **P1** | Quantum Mutual Information Reveals Energy Transfer Pathways in the FMO Complex | *J. Chem. Phys.* or *PRX* | QuTiP Lindblad on FMO | Month 6 |
| **P2** | Quantum Darwinism in Photosynthetic Energy Transfer | *New J. Phys.* or *Quantum* | Redundancy/einselection from QuTiP | Month 9 |
| **P3** | Thermodynamic Cost of Quantum Coherence in Photosynthesis | *Phys. Rev. E* or *PRX Life* | Entropy production from Lindblad/HEOM | Month 10 |
| **P4** | Machine Learning Prediction of Quantum Transport from Protein Structure | *npj Quantum Information* | GNN + QuTiP training data | Month 12 |

### Author Strategy

- Single author (Dhiraj Kumar) is sufficient for purely theoretical work.
- For experimental claims, consider a collaborator from the Gassab, Dong, or Firmenich networks.

---

## 3. Code Development Plan

### Module Architecture

```
quantum_biological_toolkit/
│
├── core/
│   ├── hamiltonian.py            # Build Trp/FMO Hamiltonians from PDB distances
│   ├── lindblad_solver.py        # QuTiP mesolve wrapper with dielectric scaling
│   ├── entropy_metrics.py        # von Neumann, QMI, Holevo, coherent info
│   └── decoherence_models.py     # Drude-Lorentz, ENAQT dephasing rates
│
├── pdb_tools/
│   ├── trp_extractor.py          # Extract Trp XYZ from PDB, compute distances
│   └── pdb_fetcher.py            # Auto-download from RCSB by PDB ID
│
├── analysis/
│   ├── channel_capacity.py       # Holevo capacity via numerical optimization
│   ├── quantum_darwinism.py      # Redundancy R(δ), pointer states, SBS
│   └── spectroscopy.py           # Simulate transient absorption signals
│
├── papers/
│   ├── paper_0_sub_tubulin/
│   ├── paper_1_qmi_fmo/
│   ├── paper_2_qd_biology/
│   ├── paper_3_thermo_cost/
│   └── paper_4_ml_transport/
│
├── quantum_optical_gateway.py        # Done — Publication-grade Lindblad
└── quantum_optical_gateway_simple.py  # Done — Pedagogical NumPy
```

### Development Phases

| Phase | Duration | Deliverables | Key Files |
|-------|----------|-------------|-----------|
| **P0** Theory | Weeks 1-2 | Sub-tubulin paper draft | `papers/p0_sub_tubulin/` |
| **P1a** PDB pipeline | Weeks 2-3 | Trp extractor, NMDA/NaV channel coordinates | `pdb_tools/trp_extractor.py` |
| **P1b** FMO baseline | Weeks 3-4 | Reproduce ENAQT curve, QMI matrix | `core/lindblad_solver.py` |
| **P2** Darwinism | Weeks 5-8 | Redundancy curves for FMO | `analysis/quantum_darwinism.py` |
| **P3** Thermodynamics | Weeks 9-10 | Entropy production rates | `analysis/channel_capacity.py` |
| **P4** ML integration | Weeks 11-16 | GNN training, feature importance | `ml/` |

---

## 4. Immediate Next Steps

### ~~Step 1 — PDB Tryptophan Extraction Script (1-2 days)~~ ✅ COMPLETE

- `src/pdb_tools/trp_extractor.py` — fetches PDB, extracts Trp CG coordinates, computes distance matrix
- Classifies pairs: "coupled" (< 1.5 nm) or "optical relay" (> 1.5 nm)
- `src/pdb_tools/batch_processor.py` — batch analysis across 63 curated membrane proteins

### ~~Step 2 — Real Hamiltonian from PDB Distances (2-3 days)~~ ✅ COMPLETE

- `QuantumOpticalGateway` updated with PDB-based Hamiltonian (J_ij ∝ 1/R_ij³)
- `HamiltonianEngine` implements H_total = H_PDB + H_clock(A)
- Clock drive adds J_clock(R) = A × 400 × exp(-R/15) cm⁻¹
- All 8 PDB targets breach classical bound (S > 2.0) under full clock drive

### ~~Step 3 — Connect Simulation to Document (1 day)~~ ✅ COMPLETE

- Section 11.5 (Geometrically Sub-Threshold Quantum Vault) added to reference doc
- `README.md` updated with full repository structure and research summary
- `RESEARCH_SUMMARY.md` written with complete results tables
- `PROJECT_PLAN.md` updated with all completed milestones
- Physical Review E manuscript at `papers/sub_tubulin_manuscript/manuscript.tex` with 4 figures
- All source modules documented with docstrings and references

### Decision Point ✅ RESOLVED — P0 first

The Sub-Tubulin theory paper (P0) is complete as a manuscript draft targeting Physical Review E. The decision was made to establish the novel hypothesis first as a pure theoretical framework before pursuing FMO-based computational papers.

## 5. Completed Milestones

| # | Milestone | Status | Deliverable |
|---|-----------|--------|-------------|
| M1 | PDB Trp extractor running on NMDA receptor | ✅ | `trp_extractor.py` + distance map |
| M2 | Hamiltonian computed from real structural data | ✅ | Updated `QuantumOpticalGateway`, `HamiltonianEngine` |
| M3 | P0 manuscript drafted | ✅ | `manuscript.tex` + 4 figures |
| M4 | Reference document complete (Parts I-XI) | ✅ | `QUANTUM_BIOLOGY_REFERENCE.md` |

## 6. Next Milestones

| # | Milestone | Timeline | Deliverable |
|---|-----------|----------|-------------|
| M5 | Submit P0 to arXiv | Next 2 weeks | Manuscript + supplementary code |
| M6 | Revise for journal submission | Month 1-2 | Address referee feedback |
| M7 | P1: QMI in FMO | Month 3-6 | Computations + manuscript |
| M8 | P2: Quantum Darwinism | Month 6-9 | Computations + manuscript |
| M9 | P3: Thermodynamic cost | Month 7-10 | Computations + manuscript |
| M10 | P4: ML for quantum transport | Month 9-12 | Computations + manuscript |
