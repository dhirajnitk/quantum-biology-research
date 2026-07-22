# Structural Constraints on Energy Transfer in Neural Tryptophan Networks

A computational analysis of the physical limits of energy transfer in neural tryptophan (Trp) networks, based on 20 PDB structures, with a two-regime characterization of non-radiative FRET vs. free-space radiative coupling.

## Preprint

A preprint of this manuscript is archived at Zenodo:
**DOI: [10.5281/zenodo.21478163](https://doi.org/10.5281/zenodo.21478163)**

## Repository Structure

```
├── papers/
│   ├── p0_sub_tubulin/               # Main manuscript (PRE target)
│   │   ├── sub_tubular_quantum.tex
│   │   └── generate_figures.py
│   ├── p1_fmo_qmi/                   # Quantum Mutual Information in FMO
│   ├── p2_quantum_darwinism/         # Quantum Darwinism in FMO
│   ├── p3_thermodynamics/            # Thermodynamic cost of coherence
│   └── p4_ml_transport/             # ML prediction of quantum transport
├── src/
│   ├── pdb_tools/                    # PDB fetching and Trp extraction
│   │   ├── trp_extractor.py
│   │   └── batch_processor.py
│   ├── core/                         # Simulation engine
│   │   ├── biophoton_relay.py        # Biophoton relay calculations
│   │   ├── quantum_optical_gateway.py
│   │   └── hamiltonian.py
│   ├── analysis/                     # Analysis modules
│   │   └── z_channel_capacity.py     # Shannon-optimal energy model
│   └── simulations/                  # Monte Carlo validation
│       ├── monte_carlo_validation.py
│       └── run_batch_analysis.py
├── DATA_SOURCES.md
└── README.md
```

## Core Result

Free-space radiative Trp-to-Trp energy transfer is infeasible (p ≈ 1.8×10⁻⁵ per exciton) due to UV tissue attenuation and spectral mismatch. Sub-1.5 nm non-radiative FRET coupling is structurally viable, with coupling strengths of J ≈ 150-300 cm⁻¹ for the closest Trp pairs in fast-signaling ion channels.

## Key Results

1. **20 PDB targets analyzed** — 18 contain 3–27 Trp residues with 1–28 quantum-coupled pairs
2. **Free-space infeasibility** — Trp-to-Trp radiative transfer p ≈ 1.8×10⁻⁵ per exciton; spectral mismatch between Trp UV emission (~320–350 nm) and acceptor NIR absorption (~600–850 nm)
3. **Two-regime analysis** — non-radiative FRET dominates at sub-1.5 nm (efficiency >0.90 for closest pairs); free-space fails beyond 2 nm
4. **Hypothetical channel constraints** — if optical coupling existed, N ≈ 5,000 cores could achieve P_success ≈ 98% with E_bit ≈ 1.85×10⁻¹⁵ J (~16× more efficient than CMOS)
5. **KCKAS contextuality scores** — static S ≤ 1.97 for all targets; classical bound S=2 not breached without active driving
6. **SDP verification** — S_max = √5 ≈ 2.236 confirmed via convex optimization

## Getting Started

```bash
# Extract Trp networks from a PDB structure
python src/pdb_tools/trp_extractor.py 7TYO --save

# Run biophoton relay calculation
python -c "from src.core.biophoton_relay import BiophotonRelay; r = BiophotonRelay(10.0); print(r.spatial_ensemble_success(5000, target_type='cco'))"

# Run Monte Carlo validation
python src/simulations/monte_carlo_validation.py

# Run Shannon-optimal energy model
python src/analysis/z_channel_capacity.py
```

## Compiling Manuscripts

```bash
cd papers/p0_sub_tubulin
pdflatex sub_tubular_quantum.tex && pdflatex sub_tubular_quantum.tex
```

## License

MIT
