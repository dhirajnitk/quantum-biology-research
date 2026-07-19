"""
trp_extractor.py  —  Extract Tryptophan networks from PDB structures.

Downloads protein structures from RCSB, locates all Trp residues,
computes the inter-Trp distance matrix, and classifies each pair
as quantum-coupled (< 1.5 nm) or optical-relay (> 1.5 nm).

Usage
-----
    python trp_extractor.py 7TYO          # NMDA receptor
    python trp_extractor.py 3ENI --chain A  # FMO complex, chain A
    python trp_extractor.py 1JFF --save    # Tubulin, save to CSV
"""

import sys
import os
import urllib.request
import numpy as np

# ── Physical threshold ──────────────────────────────────────────────
COUPLING_CUTOFF_ANGSTROM = 15.0   # 1.5 nm  —  max for Trp-Trp exciton hopping
# ─────────────────────────────────────────────────────────────────────

PDB_CACHE = os.path.join(os.path.dirname(__file__), "..", "..", "data", "cache")


def fetch_pdb(pdb_id, cache=True):
    """Download PDB file from RCSB or load from local cache."""
    os.makedirs(PDB_CACHE, exist_ok=True)
    path = os.path.join(PDB_CACHE, f"{pdb_id.upper()}.pdb")
    if os.path.exists(path):
        print(f"[+] Loaded from cache: {path}")
        return open(path).read()
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"[+] Fetching {url} ...")
    try:
        data = urllib.request.urlopen(url, timeout=30).read().decode("utf-8")
        if cache:
            with open(path, "w") as f:
                f.write(data)
        return data
    except Exception as e:
        print(f"[!] RCSB fetch failed: {e}")
        return None


def extract_trp_coordinates(pdb_text, chain_id=None):
    """Return dict: {residue_number: np.array([x, y, z])} for all Trp CG atoms."""
    centers = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        res_seq = int(line[22:26])
        chain = line[21].strip()
        if chain_id and chain != chain_id:
            continue
        # Use CG atom as the indole ring centre of mass proxy
        if res_name == "TRP" and atom_name == "CG":
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            centers[res_seq] = np.array([x, y, z])
    return centers


def distance_matrix(centers):
    """Compute Euclidean distance matrix (Å) from coordinate dict."""
    keys = sorted(centers.keys())
    n = len(keys)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            d = np.linalg.norm(centers[keys[i]] - centers[keys[j]])
            D[i, j] = D[j, i] = d
    return D, keys


def classify_pairs(D, keys):
    """Classify each Trp-Trp pair into coupling regime."""
    n = len(keys)
    coupled, relay = [], []
    for i in range(n):
        for j in range(i + 1, n):
            d = D[i, j]
            entry = (keys[i], keys[j], d)
            if d < COUPLING_CUTOFF_ANGSTROM:
                coupled.append(entry)
            else:
                relay.append(entry)
    return coupled, relay


def build_hamiltonian(D, keys):
    """Construct tight-binding Hamiltonian [cm⁻¹] from Trp distances.

    Coupling:  J_ij = J₀ * (R₀ / R_ij)³   (Dexter-like)
    Site energies:  E_i = baseline + disorder
    """
    n = len(keys)
    J0 = -80.0        # reference coupling at R0 [cm⁻¹]
    R0 = 10.0         # reference distance [Å]
    H = np.zeros((n, n))
    for i in range(n):
        E_disorder = np.random.normal(0, 30)   # static disorder
        H[i, i] = float(i * 80) + E_disorder   # site energy ramp
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] < 1e-6:
                continue
            J_ij = J0 * (R0 / D[i, j]) ** 3
            H[i, j] = H[j, i] = J_ij
    return H, keys


# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python trp_extractor.py <PDB_ID> [--chain X] [--save]")
        sys.exit(1)

    pdb_id = sys.argv[1]
    chain_id = None
    save_csv = False

    if "--chain" in sys.argv:
        idx = sys.argv.index("--chain")
        chain_id = sys.argv[idx + 1]
    if "--save" in sys.argv:
        save_csv = True

    pdb_text = fetch_pdb(pdb_id)
    if not pdb_text:
        sys.exit(1)

    centers = extract_trp_coordinates(pdb_text, chain_id)
    if not centers:
        print(f"[!] No Tryptophan residues found in {pdb_id} (chain {chain_id or 'all'}).")
        sys.exit(0)

    print(f"\n{'=' * 60}")
    print(f"  TRYPTOPHAN NETWORK ANALYSIS  -  {pdb_id.upper()}")
    print(f"  {len(centers)} Trp residues identified")
    print(f"{'=' * 60}\n")

    print("  Coordinates:")
    for seq, xyz in centers.items():
        print(f"    Trp-{seq:>4}:  ({xyz[0]:>8.3f}, {xyz[1]:>8.3f}, {xyz[2]:>8.3f})")

    D, keys = distance_matrix(centers)
    coupled, relay = classify_pairs(D, keys)

    print("\n  --- INTER-TRP DISTANCES ---")
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            d = D[i, j]
            tag = "COUPLED" if d < COUPLING_CUTOFF_ANGSTROM else "RELAY"
            print(f"    Trp-{keys[i]} -> Trp-{keys[j]}: {d:>6.2f} A  ({d*0.1:.2f} nm)  [{tag}]")

    print(f"\n  Pair summary:")
    print(f"    Quantum-coupled (< 1.5 nm):  {len(coupled)} pairs")
    print(f"    Optical relay (> 1.5 nm):    {len(relay)} pairs")

    # Build and print Hamiltonian
    H, h_keys = build_hamiltonian(D, keys)
    print(f"\n  Tight-binding Hamiltonian [cm^-1]  ({len(h_keys)}x{len(h_keys)}):")
    for row in H:
        print(f"    [{', '.join(f'{x:>8.1f}' for x in row)}]")

    if save_csv:
        import csv
        out_dir = os.path.join(os.path.dirname(__file__), "..", "..", "data")
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f"{pdb_id.upper()}_trp_network.csv")
        with open(out_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Trp_i", "Trp_j", "Distance_A", "Regime"])
            for ki, kj, d in coupled + relay:
                regime = "coupled" if d < COUPLING_CUTOFF_ANGSTROM else "relay"
                w.writerow([ki, kj, f"{d:.2f}", regime])
        print(f"\n[+] Saved pair data to {out_path}")
