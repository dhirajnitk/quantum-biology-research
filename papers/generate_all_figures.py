"""
generate_all_figures.py  —  Generate figures for P1-P4 manuscripts.

Outputs PDF files to each paper's figures/ directory.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import exp, sqrt, log, pi
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from src.core.lindblad_solver import FmoLindbladSolver

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'axes.labelsize': 11,
                     'axes.titlesize': 11, 'legend.fontsize': 9, 'figure.dpi': 300})

BLUE, ORANGE, GREEN, RED, PURPLE, GRAY = '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#7f7f7f'


def fig_p1():
    """P1: ENAQT curve + QMI heatmap."""
    out = "papers/p1_fmo_qmi/figures"
    solver = FmoLindbladSolver(disorder_std=50.0)

    # ENAQT curve
    gammas, effs = solver.enaqt_curve(np.logspace(-1, 3, 50))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.2))

    ax1.semilogx(gammas, np.array(effs) * 100, color=BLUE, lw=2.5)
    ax1.axvline(100, color=GREEN, ls='--', lw=1, alpha=0.5)
    ax1.axhline(29.4, color=GREEN, ls=':', lw=1, alpha=0.5)
    ax1.fill_between(gammas, 0, np.array(effs) * 100, alpha=0.1, color=BLUE)
    ax1.set_xlabel('Dephasing rate (cm⁻¹)')
    ax1.set_ylabel('Transfer efficiency (%)')
    ax1.set_title('(a) ENAQT with static disorder', fontweight='bold')
    ax1.set_xlim(0.1, 1000)
    ax1.grid(True, alpha=0.3)
    ax1.annotate('Anderson\nlocalization', xy=(0.2, 10), fontsize=8, color=GRAY)
    ax1.annotate('ENAQT\npeak', xy=(100, 31), fontsize=8, color=GREEN, ha='center')
    ax1.annotate('Zeno', xy=(800, 18), fontsize=8, color=RED)

    # QMI heatmap at optimal dephasing
    _, hist = solver.run_time_evolution(175 * 0.0188, t_max=10.0)
    rho = hist[-1]
    pops = np.real(np.diag(rho))
    pops = np.clip(pops, 1e-12, 1.0)
    n = 7
    qmi = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            rho_sub = np.array([[pops[i], rho[i, j]], [rho[j, i], pops[j]]])
            ev = np.clip(np.linalg.eigvalsh(rho_sub), 1e-12, 1.0)
            sj = -np.sum(ev * np.log2(ev))
            si = -pops[i] * np.log2(pops[i]) - (1-pops[i]) * np.log2(1-pops[i])
            sj_s = -pops[j] * np.log2(pops[j]) - (1-pops[j]) * np.log2(1-pops[j])
            qmi[i, j] = qmi[j, i] = max(0, si + sj_s - sj)

    im = ax2.imshow(qmi, cmap='Blues', vmin=0, vmax=0.8)
    for i in range(7):
        for j in range(7):
            ax2.text(j, i, f'{qmi[i,j]:.2f}', ha='center', va='center', fontsize=7)
    ax2.set_xticks(range(7)); ax2.set_yticks(range(7))
    ax2.set_xticklabels(range(1, 8)); ax2.set_yticklabels(range(1, 8))
    ax2.set_xlabel('Site j'); ax2.set_ylabel('Site i')
    ax2.set_title('(b) QMI matrix (bits)', fontweight='bold')
    plt.colorbar(im, ax=ax2, shrink=0.8)

    plt.tight_layout()
    plt.savefig(f'{out}/fig1_p1.pdf', bbox_inches='tight')
    plt.close()
    print(f'P1 figure saved')


def fig_p2():
    """P2: PIC curves + pointer state spectrum."""
    out = "papers/p2_quantum_darwinism/figures"
    solver = FmoLindbladSolver(disorder_std=50.0)
    _, hist = solver.run_time_evolution(175 * 0.0188, t_max=10.0)

    from src.analysis.quantum_darwinism import QuantumDarwinism
    qd = QuantumDarwinism(hist, system_site=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))

    # PIC curves
    colors = [BLUE, ORANGE, GREEN, RED]
    for t_idx, c, label in [(50, colors[0], 't=0.5 ps'), (200, colors[1], 't=2.0 ps'),
                              (500, colors[2], 't=5.0 ps'), (999, colors[3], 't=10 ps')]:
        if t_idx >= len(hist): continue
        profile = qd.compute_redundancy_curve(t_idx)
        sizes = [p[0] for p in profile]
        mis = [p[1] for p in profile]
        ax1.plot(sizes, mis, 'o-', color=c, lw=2, label=label)
    ax1.set_xlabel('Fragment size (environment sites)')
    ax1.set_ylabel('I(sys : frag) (bits)')
    ax1.set_title('(a) Partial Information Curves', fontweight='bold')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Pointer state eigenvalues
    ps, vals = qd.find_pointer_states(500)
    ax2.bar(range(1, 7), vals[:6], color=BLUE, ec='k', lw=0.5)
    ax2.set_xlabel('Eigenvalue index')
    ax2.set_ylabel('Eigenvalue')
    ax2.set_title('(b) Pointer state spectrum', fontweight='bold')
    ax2.grid(True, alpha=0.2, axis='y')

    plt.tight_layout()
    plt.savefig(f'{out}/fig2_p2.pdf', bbox_inches='tight')
    plt.close()
    print('P2 figure saved')


def fig_p3():
    """P3: Entropy trajectory + multi-dephasing comparison."""
    out = "papers/p3_thermodynamics/figures"
    solver = FmoLindbladSolver(disorder_std=50.0)
    eff, hist = solver.run_time_evolution(100 * 0.0188, t_max=10.0)

    from src.analysis.channel_capacity import FmoThermodynamicEngine
    thermo = FmoThermodynamicEngine(hist, solver.H)
    S_t = thermo.compute_entropy_trajectory()
    dS_t = thermo.compute_entropy_production_trajectory()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.5))

    t = np.arange(len(S_t)) * 0.01
    ax1.plot(t, S_t, color=BLUE, lw=2, label='S(t)')
    ax1b = ax1.twinx()
    ax1b.plot(t, dS_t, color=RED, lw=2, ls='--', label='dS/dt')
    ax1.set_xlabel('Time (ps)'); ax1.set_ylabel('Entropy S (bits)', color=BLUE)
    ax1b.set_ylabel('dS/dt (kB/ps)', color=RED)
    ax1.set_title('(a) Entropy trajectory at γ=100 cm⁻¹', fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Landauer comparison bar chart
    categories = ['Landauer\nlimit', 'Trp-CCO\nChannel', 'FMO\ncomplex', 'Classical\nCMOS']
    values = [1, 6.23e5, 3.9e2, 1e7]
    colors = [GREEN, BLUE, ORANGE, RED]
    bars = ax2.bar(categories, values, color=colors, ec='k', lw=0.5, width=0.6)
    ax2.set_yscale('log')
    ax2.set_ylabel('Energy per bit (× Landauer limit)')
    ax2.set_title('(b) Energy efficiency comparison', fontweight='bold')
    ax2.grid(True, alpha=0.2, axis='y')
    for bar, val in zip(bars, values):
        ax2.text(bar.get_x()+bar.get_width()/2, bar.get_height()*1.3,
                 f'{val:.0e}×', ha='center', fontsize=8, fontweight='bold')
    ax2.set_ylim(0.5, 2e8)

    plt.tight_layout()
    plt.savefig(f'{out}/fig3_p3.pdf', bbox_inches='tight')
    plt.close()
    print('P3 figure saved')


def fig_p4():
    """P4: RandomForest training + predictions scatter."""
    out = "papers/p4_ml_transport/figures"
    from sklearn.ensemble import RandomForestRegressor
    from src.ml.fast_transport_predictor import generate_dataset

    np.random.seed(42)
    X, y = generate_dataset(20000, seed=42)
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    rf = RandomForestRegressor(n_estimators=200, min_samples_split=2, min_samples_leaf=1,
                               random_state=42, n_jobs=-1)
    rf.fit(X_train, y_train[:, 0])
    y_pred = rf.predict(X_test)

    r2_e = 1 - np.sum((y_test[:, 0] - y_pred)**2) / max(np.sum((y_test[:, 0] - y_test[:, 0].mean())**2), 1e-10)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.5))

    # Feature importance bar chart
    importances = rf.feature_importances_
    feat_names = ['Mean coupling', 'Max coupling', 'Coupling spread', 'N sites',
                  'Dipole corr.', 'Sink distance', 'Mean distance', 'Dephasing', 'Cpl/Dep ratio']
    idx = np.argsort(importances)[::-1]
    ax1.barh(range(9), importances[idx], color=BLUE, ec='k', lw=0.5)
    ax1.set_yticks(range(9))
    ax1.set_yticklabels([feat_names[i] for i in idx], fontsize=8)
    ax1.set_xlabel('Feature importance')
    ax1.set_title('(a) RandomForest feature importance', fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='x')

    ax2.scatter(y_test[:, 0], y_pred, s=5, alpha=0.4, color=BLUE)
    ax2.plot([0, 1], [0, 1], '--', color=RED, lw=1)
    ax2.text(0.05, 0.9, f'R² = {r2_e:.3f}', fontsize=10, transform=ax2.transAxes)
    ax2.set_xlabel('True ENAQT efficiency'); ax2.set_ylabel('Predicted ENAQT')
    ax2.set_title('(b) RF predictions vs targets', fontweight='bold')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{out}/fig4_p4.pdf', bbox_inches='tight')
    plt.close()
    print('P4 figure saved')


if __name__ == "__main__":
    fig_p1(); fig_p2(); fig_p3(); fig_p4()
    print('All figures generated.')
