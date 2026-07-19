"""
manuscript_compiler.py  —  Multi-structure comparative analysis for publication.

Integrates all analysis modules (Trp extraction, spatial ensemble, KCKAS
contextuality, CHSH Bell test, Z-channel capacity) into a single
manuscript-ready data matrix.

Outputs the tables and figures needed for the Results section of a
paper targeting Physical Review E or J. Chem. Phys.

References
----------
[B2024] Babcock et al. JPCB (2024)
[KCKAS2008] Klyachko et al. PRL 101 (2008)
[CHSH1969] Clauser et al. PRL 23 (1969)
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.pdb_tools.trp_extractor import fetch_pdb, extract_trp_coordinates, distance_matrix
from src.core.biophoton_relay import BiophotonRelay, channel_capacity, THERMAL_SUPPRESSION
from src.analysis.pdb_contextuality import compute_kckas_from_pdb, KCKAS_CLASSICAL_BOUND, KCKAS_QUANTUM_BOUND
from src.analysis.pdb_bell_test import compute_chsh_from_pdb
from src.analysis.z_channel_capacity import ensemble_gap_to_shannon

# ── Master target list ─────────────────────────────────────────────
# Each entry: PDB ID → metadata for the manuscript table
TARGETS = [
    # (pdb_id, name, class, dielectric)
    ("6CNO",  "NMDA receptor GluN1/GluN2B",  "Ligand-gated ion channel", 2.1),
    ("6LQA",  "NaV1.4 voltage-gated sodium",  "Voltage-gated ion channel", 2.3),
    ("7KOX",  "Alpha-7 nicotinic AChR",       "Ligand-gated ion channel", 2.0),
    ("7TYO",  "NMDA receptor GluN1/GluN2A",  "Ligand-gated ion channel", 2.1),
    ("6J8J",  "NaV1.4 (alternative)",         "Voltage-gated ion channel", 2.3),
    ("6PV7",  "Nicotinic ACh receptor",       "Ligand-gated ion channel", 2.0),
    ("1BL8",  "KcsA potassium channel",       "Voltage-gated ion channel", 2.0),
    ("1YAG",  "SNARE complex",                "Synaptic vesicle anchor",  2.6),
    ("1JFF",  "Tubulin dimer (bovine)",        "Cytoskeletal filament",    2.1),
    ("3N2K",  "Tubulin dimer (mammalian)",     "Cytoskeletal filament",    2.1),
]


def analyse_target(pdb_id, name, pdb_class, dielectric):
    """Run full analysis pipeline on one PDB target.

    Returns a dict of all metrics for the manuscript table.
    """
    text = fetch_pdb(pdb_id)
    if not text:
        return None

    centres = extract_trp_coordinates(text)
    n_trp = len(centres)

    # Trp pair analysis
    coupled_pairs = 0
    relay_pairs = 0
    mean_distance = 0.0
    if n_trp >= 2:
        D, keys = distance_matrix(centres)
        n = len(keys)
        coupled_pairs = sum(1 for i in range(n) for j in range(i+1, n) if D[i, j] < 15.0)
        relay_pairs = n * (n - 1) // 2 - coupled_pairs
        mean_distance = float(np.mean([D[i, j] for i in range(n) for j in range(i+1, n)]))

    # KCKAS contextuality
    kckas = compute_kckas_from_pdb(pdb_id)
    kckas_S = kckas["kckas_sum"] if kckas else 0.0
    kckas_F = kckas["coherence_factor"] if kckas else 0.0

    # CHSH Bell test
    chsh = compute_chsh_from_pdb(pdb_id)
    chsh_S = chsh["S_value"] if chsh else 0.0

    # Spatial ensemble (CCO target)
    relay = BiophotonRelay(mean_distance * 0.1 if mean_distance > 0 else 3.0)
    ensemble_sizes = [1000, 5000, 10000]
    results, p_per_core = relay.compute_ensemble_scan(ensemble_sizes, target_type="cco")
    p_ensemble = {r["n_cores"]: r["p_ensemble"] for r in results}

    # Critical N for 95% power
    N_95 = int(np.ceil(np.log(0.05) / np.log(1.0 - p_per_core))) if p_per_core > 0 else 0

    # Z-channel capacity at N=5,000
    cap_info = ensemble_gap_to_shannon(p_per_core, 5000)

    return {
        "pdb_id": pdb_id,
        "name": name,
        "class": pdb_class,
        "n_trp": n_trp,
        "n_tyr_phe": 0,  # Would need full PDB parsing for Tyr/Phe counts
        "coupled_pairs": coupled_pairs,
        "relay_pairs": relay_pairs,
        "mean_distance_A": mean_distance,
        "dielectric": dielectric,
        "p_per_core": p_per_core,
        "p_1000": p_ensemble.get(1000, 0),
        "p_5000": p_ensemble.get(5000, 0),
        "p_10000": p_ensemble.get(10000, 0),
        "N_95": N_95,
        "kckas_S": kckas_S,
        "kckas_F_coh": kckas_F,
        "kckas_violation": kckas_S - KCKAS_CLASSICAL_BOUND if kckas_S > KCKAS_CLASSICAL_BOUND else 0.0,
        "chsh_S": chsh_S,
        "z_capacity_5000": cap_info["c_ensemble"] if cap_info else 0.0,
        "z_gap_to_shannon": cap_info["gap"] if cap_info else 0.0,
    }


def compile_manuscript_table(targets):
    """Analyse all targets and format the manuscript data matrix."""
    print("=" * 100)
    print("  MANUSCRIPT DATA MATRIX — Comparative Structural Analysis")
    print("  Targets: ligand-gated ion channels, voltage-gated, synaptic, cytoskeletal")
    print("=" * 100)

    # ── Table 1: Structural properties ──────────────────────────
    print(f"\n  TABLE 1: Structural Properties of Neural Trp Networks")
    print(f"  {'PDB':<6} {'Name':<30} {'Class':<30} {'Trp':<5} {'Coupled':<8} {'Relay':<8} {'Dist(A)':<8} {'Dielec':<8}")
    print(f"  {'-'*6} {'-'*30} {'-'*30} {'-'*5} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")

    results = []
    for pdb_id, name, pdb_class, eps in targets:
        r = analyse_target(pdb_id, name, pdb_class, eps)
        if r:
            results.append(r)
            print(f"  {r['pdb_id']:<6} {r['name'][:28]:<30} {r['class'][:28]:<30} {r['n_trp']:<5} {r['coupled_pairs']:<8} {r['relay_pairs']:<8} {r['mean_distance_A']:<8.1f} {r['dielectric']:<8.1f}")

    # ── Table 2: Quantum signatures ─────────────────────────────
    print(f"\n  TABLE 2: Quantum Contextuality & Bell-CHSH Signatures")
    print(f"  {'PDB':<6} {'KCKAS S':<10} {'F_coh':<8} {'Violation':<12} {'CHSH S':<10} {'Status':<20}")
    print(f"  {'-'*6} {'-'*10} {'-'*8} {'-'*12} {'-'*10} {'-'*20}")

    for r in results:
        kckas_status = "CONTEXTUAL" if r['kckas_S'] > KCKAS_CLASSICAL_BOUND else "Classical-like"
        chsh_status = "Bell-violating" if r['chsh_S'] > 2.0 else "Classical"
        status = f"{kckas_status}, {chsh_status}"
        print(f"  {r['pdb_id']:<6} {r['kckas_S']:<10.4f} {r['kckas_F_coh']:<8.4f} {r['kckas_violation']:<12.4f} {r['chsh_S']:<10.4f} {status:<20}")

    # ── Table 3: Spatial ensemble gating ─────────────────────────
    print(f"\n  TABLE 3: Z-Channel Gating — Spatial Ensemble Thresholds")
    print(f"  {'PDB':<6} {'p/core':<12} {'N_95%':<10} {'P@1k':<10} {'P@5k':<10} {'P@10k':<12} {'C(5k)':<10} {'Gap to Sh.':<12}")
    print(f"  {'-'*6} {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*12} {'-'*10} {'-'*12}")

    for r in results:
        print(f"  {r['pdb_id']:<6} {r['p_per_core']:<12.2e} {r['N_95']:<10} {r['p_1000']*100:<8.1f}% {r['p_5000']*100:<8.1f}% {r['p_10000']*100:<8.1f}% {r['z_capacity_5000']:<10.4f} {r['z_gap_to_shannon']:<12.6f}")

    # ── Summary statistics ──────────────────────────────────────
    print(f"\n  SUMMARY")
    n_contextual = sum(1 for r in results if r['kckas_S'] > KCKAS_CLASSICAL_BOUND)
    n_total = len(results)
    mean_coupled = float(np.mean([r['coupled_pairs'] for r in results]))
    mean_N95 = float(np.mean([r['N_95'] for r in results]))
    print(f"    Targets analysed:     {n_total}")
    print(f"    Mean coupled pairs:   {mean_coupled:.1f}")
    print(f"    Mean N for 95% power: {mean_N95:.0f}")
    print(f"    Contextual (S>2):     {n_contextual}/{n_total}")

    print(f"\n  NOTE: KCKAS bound S_classical = {KCKAS_CLASSICAL_BOUND:.0f}, S_quantum_max = {KCKAS_QUANTUM_BOUND:.4f}")
    print(f"        CHSH bound S_classical = 2.0, S_quantum_max = 2.828")
    print(f"        All optical params from [B2024], dielectrics from literature.")
    print("=" * 100)

    return results


if __name__ == "__main__":
    results = compile_manuscript_table(TARGETS)
