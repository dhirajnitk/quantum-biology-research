"""
generate_figures.py  —  Generate all figures for the PRE manuscript.

Outputs PDF files to papers/sub_tubulin_manuscript/figures/

Requires: numpy, matplotlib
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import exp, sqrt, log, pi
import os

OUT_DIR = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Style ───────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 300,
})

BLUE = '#1f77b4'
ORANGE = '#ff7f0e'
GREEN = '#2ca02c'
RED = '#d62728'
PURPLE = '#9467bd'
GRAY = '#7f7f7f'


# ═══════════════════════════════════════════════════════════════════
#  Figure 1: Trp-CCO Channel diagram + spatial ensemble summation
# ═══════════════════════════════════════════════════════════════════

def fig1_trp_cco_channel():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.5, 3.2))

    # Panel A: Trp-CCO Channel transition diagram (conceptual)
    ax1.set_xlim(-0.3, 2.3)
    ax1.set_ylim(-0.3, 1.8)
    ax1.set_aspect('equal')
    ax1.axis('off')

    # Nodes
    circle0 = plt.Circle((0, 1.0), 0.25, color=BLUE, ec='k', lw=2)
    circle1 = plt.Circle((2, 1.0), 0.25, color=ORANGE, ec='k', lw=2)
    circle0b = plt.Circle((0, 0), 0.25, color=BLUE, ec='k', lw=2)
    circle1b = plt.Circle((2, 0), 0.25, color=ORANGE, ec='k', lw=2, alpha=0.3)
    ax1.add_patch(circle0)
    ax1.add_patch(circle1)
    ax1.add_patch(circle0b)
    ax1.add_patch(circle1b)

    ax1.text(0, 1.0, '0', ha='center', va='center', fontsize=14, fontweight='bold', color='white')
    ax1.text(2, 1.0, '1', ha='center', va='center', fontsize=14, fontweight='bold', color='white')
    ax1.text(0, 0, '0', ha='center', va='center', fontsize=14, fontweight='bold', color='white')
    ax1.text(2, 0, '1', ha='center', va='center', fontsize=14, fontweight='bold', color='gray', alpha=0.3)

    # Transitions
    ax1.annotate('', xy=(0.25, 1.0), xytext=(0, 1.0),
                 arrowprops=dict(arrowstyle='->', lw=2, color=GRAY))
    ax1.text(0.1, 1.08, 'p=1', fontsize=9, color=GRAY)

    ax1.annotate('', xy=(1.75, 0.75), xytext=(2.25, 1.0),
                 arrowprops=dict(arrowstyle='->', lw=1.5, color=GRAY, connectionstyle='arc3,rad=-0.5'))
    ax1.text(2.1, 0.6, 'p≈0.06%', fontsize=8, color=GRAY, ha='center')

    ax1.annotate('', xy=(1.75, 0.05), xytext=(2.25, 1.0),
                 arrowprops=dict(arrowstyle='->', lw=1.5, color=GRAY, connectionstyle='arc3,rad=0.5', linestyle='dashed'))
    ax1.text(2.1, 0.15, '1-p≈99.94%', fontsize=8, color=GRAY, ha='center')

    ax1.text(-0.1, 1.35, 'Input', fontsize=10, fontweight='bold')
    ax1.text(1.9, 1.35, 'Output', fontsize=10, fontweight='bold')
    ax1.set_title('(a) Trp-CCO Channel', fontsize=11, fontweight='bold')

    # Panel B: Spatial ensemble P(N)
    p = 7.84e-4
    N = np.logspace(0, 5, 200).astype(int)
    P = 1 - (1 - p) ** N

    ax2.semilogx(N, P * 100, color=BLUE, lw=2.5)
    ax2.axhline(95, color=GREEN, ls='--', lw=1, alpha=0.7, label='95% threshold')
    ax2.axhline(99, color=RED, ls='--', lw=1, alpha=0.7, label='99% threshold')
    ax2.axvline(5060, color=GREEN, ls=':', lw=1, alpha=0.5)
    ax2.axvline(7780, color=RED, ls=':', lw=1, alpha=0.5)

    ax2.fill_between(N, 0, P * 100, alpha=0.1, color=BLUE)
    ax2.set_xlabel('Spatial ensemble size N (cores)')
    ax2.set_ylabel('Gating success probability (%)')
    ax2.set_title('(b) Spatial ensemble summation', fontsize=11, fontweight='bold')
    ax2.set_xlim(1, 1e5)
    ax2.set_ylim(0, 100)
    ax2.legend(loc='lower right', framealpha=0.9)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig1_trp_cco_channel.pdf')
    plt.savefig(path, bbox_inches='tight')
    print(f"Saved {path}")
    plt.close()


# ═══════════════════════════════════════════════════════════════════
#  Figure 2: Geometry tax — KCKAS S vs static coherence
# ═══════════════════════════════════════════════════════════════════

def fig2_geometry_tax():
    fig, ax = plt.subplots(figsize=(6, 4))

    targets = ['1BL8\n(KcsA)', '6PV7\n(nAChR)', '7TYO\n(NMDA)', '6LQA\n(NaV)',
               '6CNO\n(NMDA)', '7KOX\n(nAChR)', '6J8J\n(NaV)', '1YAG\n(SNARE)']
    S_static = [1.9735, 1.9174, 1.8092, 1.7236, 1.3504, 1.5824, 1.7323, 1.3124]
    F_coh = [0.9868, 0.9587, 0.9046, 0.8618, 0.6752, 0.7912, 0.8662, 0.6562]
    S_clock1 = [2.2274, 2.1658, 2.1536, 2.1051, 2.0987, 2.1156, 2.0616, 2.1299]
    A_crit = [0.03, 0.20, 0.36, 0.57, 0.71, 0.61, 0.71, 0.63]

    x = np.arange(len(targets))
    width = 0.28

    bars1 = ax.bar(x - width, S_static, width, color=BLUE, ec='k', lw=0.5, label='Static geometry')
    bars2 = ax.bar(x, S_clock1, width, color=ORANGE, ec='k', lw=0.5, label='With clock drive (A=1)')

    ax.axhline(2.0, color=RED, ls='--', lw=2, alpha=0.8, label='Classical bound (S=2)')
    ax.axhline(2.2361, color=PURPLE, ls=':', lw=1.5, alpha=0.6, label='Quantum maximum (√5)')

    # Annotate A_crit above bars
    for i, (a, s) in enumerate(zip(A_crit, S_clock1)):
        ax.text(i, s + 0.05, f'Ac={a:.2f}', ha='center', va='bottom', fontsize=7, color=GRAY, rotation=45)

    ax.set_xticks(x)
    ax.set_xticklabels(targets, fontsize=8)
    ax.set_ylabel('KCKAS contextuality sum S')
    ax.set_title('Geometry tax: symmetric channels need less clock drive', fontsize=11, fontweight='bold')
    ax.legend(loc='lower right', framealpha=0.9, fontsize=8)
    ax.set_ylim(0, 2.5)
    ax.grid(True, alpha=0.2, axis='y')

    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig2_geometry_tax.pdf')
    plt.savefig(path, bbox_inches='tight')
    print(f"Saved {path}")
    plt.close()


# ═══════════════════════════════════════════════════════════════════
#  Figure 3: Hamiltonian engine — S(A) curves
# ═══════════════════════════════════════════════════════════════════

def fig3_hamiltonian():
    fig, ax = plt.subplots(figsize=(6, 4))

    pdb_data = {
        '1BL8 (KcsA)':  {'S0': 1.9735, 'F': 0.9868, 'A_crit': 0.03},
        '6PV7 (nAChR)': {'S0': 1.9174, 'F': 0.9587, 'A_crit': 0.20},
        '7TYO (NMDA)':  {'S0': 1.8092, 'F': 0.9046, 'A_crit': 0.36},
        '6LQA (NaV)':   {'S0': 1.7236, 'F': 0.8618, 'A_crit': 0.57},
        '6CNO (NMDA)':  {'S0': 1.3504, 'F': 0.6752, 'A_crit': 0.71},
        '7KOX (nAChR)': {'S0': 1.5824, 'F': 0.7912, 'A_crit': 0.61},
    }

    colors = [BLUE, ORANGE, GREEN, RED, PURPLE, GRAY]
    S_q = sqrt(5)
    delta_E = 80.0
    J_max = 400.0

    A = np.linspace(0, 1, 200)

    for (label, data), color in zip(pdb_data.items(), colors):
        S0 = data['S0']
        F = data['F']
        A_c = data['A_crit']

        # Simulate clock coupling: J_eff(A) ~ A * J_max * F
        J_eff = A * J_max * F
        enhancement = 1 - exp(-J_eff / delta_E)
        S = S0 + (S_q - S0) * enhancement

        ax.plot(A, S, color=color, lw=2, label=f'{label} (Ac={A_c:.2f})')
        ax.scatter([A_c], [2.0001], color=color, s=40, zorder=5, edgecolor='k')

    ax.axhline(2.0, color=RED, ls='--', lw=2, alpha=0.8)
    ax.axhline(2.2361, color=PURPLE, ls=':', lw=1.5, alpha=0.5)

    ax.annotate('Classical bound S=2', xy=(0.02, 2.02), fontsize=9, color=RED)
    ax.annotate('Quantum max √5 ≈ 2.236', xy=(0.02, 2.25), fontsize=9, color=PURPLE)

    ax.set_xlabel('Clock amplitude A')
    ax.set_ylabel('KCKAS contextuality sum S')
    ax.set_title('Clock-driven coherence rescue across neural targets', fontsize=11, fontweight='bold')
    ax.legend(loc='lower right', framealpha=0.9, fontsize=7, ncol=1)
    ax.set_ylim(1.0, 2.5)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig3_hamiltonian.pdf')
    plt.savefig(path, bbox_inches='tight')
    print(f"Saved {path}")
    plt.close()


# ═══════════════════════════════════════════════════════════════════
#  Figure 4: Energy efficiency comparison
# ═══════════════════════════════════════════════════════════════════

def fig4_energy():
    fig, ax = plt.subplots(figsize=(5, 3.5))

    categories = ['Landauer\nlimit', 'Trp-CCO\n(UV)', 'Flavin\n(Blue)', 'CCO\n(NIR)', 'FMO\ncomplex', 'CMOS']
    values = [1, 6.23e5, 3.88e5, 2.05e5, 3.9e2, 1e7]
    colors_vals = [GREEN, BLUE, '#4287f5', '#6b2d82', ORANGE, RED]

    bars = ax.bar(categories, values, color=colors_vals, ec='k', lw=1, width=0.6)

    ax.set_yscale('log')
    ax.set_ylabel('Energy per bit (× Landauer limit)')
    ax.set_title('Energy efficiency comparison', fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.2, axis='y')

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.3,
                f'{val:.0e}×', ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax.set_ylim(0.5, 2e8)

    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'fig4_energy.pdf')
    plt.savefig(path, bbox_inches='tight')
    print(f"Saved {path}")
    plt.close()


if __name__ == "__main__":
    fig1_trp_cco_channel()
    fig2_geometry_tax()
    fig3_hamiltonian()
    fig4_energy()
    print("\nAll figures generated in", OUT_DIR)
