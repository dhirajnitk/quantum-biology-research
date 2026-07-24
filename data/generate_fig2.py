"""Regenerate Figure 2 for the manuscript with real PDB data."""
import sys, os, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, 'src')
from pdb_tools.trp_extractor import fetch_pdb, extract_trp_data
from analysis.quantum_hamiltonian_engine import build_pdb_hamiltonian

pdb_ids = ['1BL8','2BG9','6IRA','6AGF','4PE5','7KOX','6PM6','5C1M','3J5P','3GD8',
           '1SFC','3KG2','2A79','5GJV','5SVK','1U19','1JFF','3N2K','6HUO']

# Collect all Trp-Trp pair distances from all 19 PDBs
all_distances = []
for pid in pdb_ids:
    text = fetch_pdb(pid, cache=True)
    if not text:
        continue
    centroids, dipoles = extract_trp_data(text)
    keys = sorted(centroids.keys())
    n = len(keys)
    for i in range(n):
        for j in range(i + 1, n):
            ri = centroids[keys[i]]
            rj = centroids[keys[j]]
            d = np.linalg.norm(np.array(ri) - np.array(rj))
            all_distances.append(d * 0.1)  # convert Angstrom to nm

all_distances = np.array(all_distances)
print(f'Total pairs collected: {len(all_distances)}')

# ── Two-panel figure ────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

# Panel (a): Two-regime efficiency
R0 = 1.3     # nm
Phi_F = 0.21
sigma = 1.8e-22  # m^2
sigma_nm2 = sigma * 1e18  # convert to nm^2
alpha = 1000    # m^-1
alpha_nm = alpha * 1e-9  # convert to nm^-1

r = np.logspace(-1, 2, 500)  # 0.1 to 100 nm
P_FRET = 1 / (1 + (r / R0)**6)
P_rad = Phi_F * sigma_nm2 / (4 * np.pi * r**2) * np.exp(-alpha_nm * r)

ax1.loglog(r, P_FRET, 'b-', linewidth=2, label='FRET')
ax1.loglog(r, P_rad, 'orange', linewidth=2, label='Free-space radiative')
ax1.axvline(1.5, color='green', linestyle='--', alpha=0.7, label='1.5 nm cutoff')
ax1.set_xlabel('Separation distance (nm)')
ax1.set_ylabel('Transfer efficiency')
ax1.set_xlim(0.1, 100)
ax1.set_ylim(1e-10, 2)
ax1.legend(fontsize=8)
ax1.set_title('(a) Two-regime energy transfer', fontsize=9)
ax1.grid(True, which='both', alpha=0.3)

# Panel (b): Pair-distance histogram
ax2.hist(all_distances, bins=50, range=(0, 10), color='steelblue', edgecolor='white', alpha=0.8)
ax2.axvline(1.5, color='green', linestyle='--', alpha=0.7, label='1.5 nm cutoff')
ax2.axvline(5.0, color='red', linestyle=':', alpha=0.5, label='5.0 nm cutoff')
ax2.set_xlabel('Trp-Trp distance (nm)')
ax2.set_ylabel('Count')
ax2.set_title('(b) Trp pair distances (19 PDB targets)', fontsize=9)
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
out_path = os.path.join(os.path.dirname(__file__), '..', 'papers', 'p0_sub_tubulin', 'figures', 'fig2_regimes.pdf')
plt.savefig(out_path, dpi=300)
print(f'Figure saved to {out_path}')
