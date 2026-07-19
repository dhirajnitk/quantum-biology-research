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
| 11.5 | Simulation Pipeline | 🔴 Needs writing |
| 11.6 | Paper Roadmap | 🔴 Needs writing |

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

### Step 1 — PDB Tryptophan Extraction Script (1-2 days)

- Use Biopython to fetch real PDB structures (target: 7TYO = NMDA receptor)
- Extract all Trp residue coordinates
- Compute inter-Trp distance matrix
- Classify each pair as "coupled" (< 1.5 nm) or "optical relay" (> 1.5 nm)

### Step 2 — Real Hamiltonian from PDB Distances (2-3 days)

- Replace the mock 4×4 Hamiltonian in `QuantumOpticalGateway` with one computed from actual Trp distances
- Coupling strength: `J_ij ∝ 1 / R_ij³` (Dexter) + dipole orientation factor

### Step 3 — Connect Simulation to Document (1 day)

- Add Section 11.5 documenting both simulation modules
- Add Section 11.6 with the paper roadmap

### Decision Point

The #1 choice is whether to write the **Sub-Tubulin theory paper (P0)** first as a pure theoretical framework, or jump straight into **FMO-based computational papers (P1-P4)** using established protein data.

- **P0 first** — Establishes your novel hypothesis; positions you as a theorist.
- **P1-P4 first** — Builds publication record with tractable simulations; defers high-risk theory.

---

## 5. Key Milestones

| # | Milestone | Date | Deliverable |
|---|-----------|------|-------------|
| M1 | PDB Trp extractor running on NMDA receptor | Week 2 | `trp_extractor.py` + distance map |
| M2 | Hamiltonian computed from real structural data | Week 3 | Updated `QuantumOpticalGateway` |
| M3 | P0 submitted to arXiv | Week 8 | Manuscript + figures |
| M4 | FMO ENAQT curve reproduced | Week 4 | Jupyter notebook |
| M5 | P1 submitted to arXiv | Month 6 | Manuscript |
| M6 | P2 submitted to arXiv | Month 9 | Manuscript |
| M7 | P3 submitted to arXiv | Month 10 | Manuscript |
| M8 | P4 submitted to arXiv | Month 12 | Manuscript + code release |
