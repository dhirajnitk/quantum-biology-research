# Sub-Tubulin Quantum Information Processing

A theoretical and computational framework for understanding how quantum-coherent information processing occurs within intramolecular Tryptophan (Trp) networks buried inside individual synaptic proteins — the **sub-tubulin hypothesis**.

## Repository Structure

```
├── docs/
│   ├── reference/QUANTUM_BIOLOGY_REFERENCE.md   # Main reference (Parts I-XI)
│   ├── architecture/                            # Narrative essay
│   └── research_plan/                           # 12-section research plan
├── src/
│   ├── pdb_tools/                               # PDB fetching and Trp extraction
│   │   ├── trp_extractor.py
│   │   └── batch_processor.py
│   └── core/                                    # Simulation engine
│       ├── quantum_optical_gateway.py           # QuTiP Lindblad simulation
│       └── quantum_optical_gateway_simple.py    # Pedagogical NumPy version
├── papers/
│   └── p0_sub_tubulin/outline.md                # P0 paper outline
├── data/
│   ├── cache/                                   # PDB structure cache (gitignored)
│   └── *_batch.csv                              # Batch analysis results
├── PROJECT_PLAN.md
└── README.md
```

## Core Hypothesis

The human brain's 20-watt efficiency cannot be explained by classical action potentials alone. Quantum-coherent information processing occurs within **1-5 nm hydrophobic protein pockets** containing dense Tryptophan aromatic networks. These serve as ultra-fast local quantum filters, communicating via near-infrared biophoton emission guided by lipid membrane waveguides — all synchronized by macro-scale electrical rhythms.

## Key Claims

1. **Trp networks** in hydrophobic protein cores maintain 50-100 fs coherence (shielded by ε ≈ 2 dielectric)
2. **Biophoton emission** (600-1300 nm) relays quantum outputs between cores via lipid waveguide
3. **Multi-scale clock hierarchy** (Hz → ps → fs) phase-synchronizes billions of independent filters
4. **Observability crossover** — Trp networks are spectroscopically verifiable (unlike nuclear spin models)

## Related Research Groups & Published Data

Several labs are producing directly relevant experimental data:

| Group | Key Finding | Year |
|-------|-------------|------|
| **Babcock, Kurian et al.** (NIH/Howard) | UV superradiance from mega-networks of Trp in microtubules & amyloids; collective quantum states survive at room temperature | 2024 |
| **Gassab, Craddock et al.** | HEOM on tubulin Trp networks (~13 fs dephasing at 310K); quantum information flow in Trp networks | 2026 |
| **Firmenich et al.** | Rigorous non-perturbative HEOM on 1JFF tryptophan network; introduced equilibrium-to-functionality gap | 2026 |
| **Patwa, Kurian et al.** | Quantum-enhanced photoprotection in Trp neuroprotein architectures; subradiant dark states | 2024 |
| **Xin et al.** | Molecular dynamics-derived coloured noise mediates ENAQT of Trp excitons in tubulin | 2026 |
| **Echternach** | Geometric foundations for collective quantum phenomena in microtubules; nanophotonic analogs | 2025-2026 |

**Key papers:**
- Babcock et al. (2024) *J. Phys. Chem. Lett.* — UV superradiance from Trp mega-networks
- Patwa et al. (2024) *Frontiers in Physics* — Quantum-enhanced photoprotection in Trp architectures
- Gassab et al. (2026) *Entropy* — Quantum information flow in microtubule Trp networks
- Firmenich et al. (2026) *bioRxiv* — HEOM on tubulin Trp networks
- Gassab & Craddock (2026) *J. Phys. Photonics* — Optical properties of cytoskeleton

## Getting Started

```bash
# Extract Trp networks from a PDB structure
python src/pdb_tools/trp_extractor.py 7TYO --save

# Batch analyse 63 curated membrane proteins
python src/pdb_tools/batch_processor.py

# Run the Quantum-Optical Gateway simulation
python src/core/quantum_optical_gateway.py
```

## Dataset

Currently analysing 63 curated membrane proteins from PDB (ion channels, GPCRs, transporters, aquaporins). The curated list spans:

- **Voltage-gated ion channels:** NaV, CaV, KV, KcsA
- **Ligand-gated ion channels:** NMDA, nAChR, GABA-A, GlyR, P2X
- **GPCRs:** Rhodopsin, β2-adrenergic, adenosine, opioid, cannabinoid
- **Transporters:** SERT, DAT, GLUT
- **Water channels:** Aquaporins
- **Junction proteins:** Connexins, TRPV

## Status

Active theoretical research. All code runs on a standard laptop with Python + QuTiP. No wet lab required.

## License

MIT
