"""
trp_extractor.py  —  Extract Tryptophan networks from PDB structures.

Downloads PDB structures from RCSB, locates all Trp residues, extracts
both centroids and transition dipole vectors from full indole ring atoms,
computes inter-Trp distances and dipole orientation factors.

The dipole-dipole coupling including orientation is:
    J_ij = J0 * (R0 / R_ij)^3 * kappa_ij
where kappa = mu_i_hat . mu_j_hat - 3 (mu_i_hat . R_hat)(mu_j_hat . R_hat)

Usage
-----
    python trp_extractor.py 7TYO
    python trp_extractor.py 3ENI --chain A --save
"""

import sys, os, urllib.request
import numpy as np

COUPLING_CUTOFF_ANGSTROM = 15.0   # 1.5 nm
PDB_CACHE = os.path.join(os.path.dirname(__file__), "..", "..", "data", "cache")

# Atoms forming the Trp indole ring (for dipole orientation)
INDOLE_ATOMS = {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}


def fetch_pdb(pdb_id, cache=True):
    """Download PDB file from RCSB or load from local cache."""
    os.makedirs(PDB_CACHE, exist_ok=True)
    path = os.path.join(PDB_CACHE, f"{pdb_id.upper()}.pdb")
    if os.path.exists(path):
        return open(path).read()
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        data = urllib.request.urlopen(url, timeout=30).read().decode("utf-8")
        if cache:
            with open(path, "w") as f:
                f.write(data)
        return data
    except Exception as e:
        print(f"[!] RCSB fetch failed for {pdb_id}: {e}")
        return None


def extract_trp_data(pdb_text, chain_id=None):
    """Extract Trp centroids AND transition dipole vectors.

    Returns
    -------
    centroids : dict {res_num: np.array([x,y,z])}
        CG atom coordinates (ring centre-of-mass proxy).
    dipoles : dict {res_num: np.array([dx,dy,dz])}
        Unit vector along the indole ring long axis (CG->CH2 direction).
    """
    ring_atoms = {}  # {res_num: [x,y,z, atom_name]}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21].strip()
        if chain_id and chain != chain_id:
            continue
        if res_name != "TRP":
            continue
        try:
            res_seq = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue

        if atom_name not in INDOLE_ATOMS:
            continue
        if res_seq not in ring_atoms:
            ring_atoms[res_seq] = []
        ring_atoms[res_seq].append((atom_name, np.array([x, y, z])))

    centroids = {}
    dipoles = {}
    for res_seq, atoms in ring_atoms.items():
        coord_dict = {name: xyz for name, xyz in atoms}
        # Centroid = CG atom
        if "CG" in coord_dict:
            centroids[res_seq] = coord_dict["CG"]
        else:
            # Fall back to mean of all ring atoms
            centroids[res_seq] = np.mean([xyz for _, xyz in atoms], axis=0)

        # Transition dipole: unit vector along the ring long axis
        # The Trp indole ring has approximate D2h symmetry.
        # The transition dipole for the S0->S1 (La) transition lies
        # along the long axis of the ring (CG->CH2 direction).
        if "CG" in coord_dict and "CH2" in coord_dict:
            d = coord_dict["CH2"] - coord_dict["CG"]
            norm = np.linalg.norm(d)
            dipoles[res_seq] = d / norm if norm > 0 else np.array([1.0, 0.0, 0.0])
        elif "CD1" in coord_dict and "CD2" in coord_dict:
            # Alternative: use the line connecting CD1 and CD2
            d = coord_dict["CD2"] - coord_dict["CD1"]
            norm = np.linalg.norm(d)
            dipoles[res_seq] = d / norm if norm > 0 else np.array([1.0, 0.0, 0.0])
        else:
            dipoles[res_seq] = np.array([1.0, 0.0, 0.0])

    return centroids, dipoles


def extract_trp_coordinates(pdb_text, chain_id=None):
    """Backward-compatible wrapper: returns only centroid dict."""
    centroids, _ = extract_trp_data(pdb_text, chain_id)
    return centroids


def distance_matrix(centers):
    """Euclidean distance matrix from centroid dict."""
    keys = sorted(centers.keys())
    n = len(keys)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            D[i, j] = D[j, i] = np.linalg.norm(centers[keys[i]] - centers[keys[j]])
    return D, keys


def dipole_orientation_factor(mu_i, mu_j, R_ij):
    """Kappa factor for dipole-dipole coupling.

    kappa = mu_i_hat . mu_j_hat - 3 (mu_i_hat . R_hat)(mu_j_hat . R_hat)

    |kappa| ranges from 0 to 2. kappa = 0 when dipoles are orthogonal
    or at the magic angle. |kappa| = 2 for collinear head-to-tail.
    """
    r_hat = R_ij / (np.linalg.norm(R_ij) + 1e-30)
    mu_i_u = mu_i / (np.linalg.norm(mu_i) + 1e-30)
    mu_j_u = mu_j / (np.linalg.norm(mu_j) + 1e-30)

    kappa = np.dot(mu_i_u, mu_j_u) - 3 * np.dot(mu_i_u, r_hat) * np.dot(mu_j_u, r_hat)
    return kappa


def classify_pairs(D, keys):
    """Classify Trp pairs by distance regime."""
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


def build_hamiltonian(D, keys, dipoles=None, dielectric=2.0, centroids=None):
    """Tight-binding Hamiltonian with dipole orientation.

    J_ij = J0 * (R0 / R_ij)^3 * kappa_ij / sqrt(dielectric)
    where kappa uses the REAL inter-residue displacement vector.

    Without dipoles: J_ij = J0 * (R0 / R_ij)^3 / sqrt(dielectric)
    """
    n = len(keys)
    J0 = -80.0     # cm^-1 at R0
    R0 = 10.0      # Angstroms
    rng = np.random.RandomState(42)

    H = np.zeros((n, n))
    for i in range(n):
        H[i, i] = float(i * 80) + rng.normal(0, 30)

    key_list = list(keys) if isinstance(keys, (list, tuple)) else keys
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] < 1e-6:
                continue
            J_ij = J0 * (R0 / D[i, j]) ** 3 / np.sqrt(dielectric)

            # Apply dipole orientation factor with REAL displacement vector
            if dipoles is not None:
                ki, kj = key_list[i], key_list[j]
                if ki in dipoles and kj in dipoles:
                    # Compute actual displacement vector if centroids available
                    if centroids is not None and ki in centroids and kj in centroids:
                        R_vec = centroids[kj] - centroids[ki]
                    else:
                        R_vec = np.array([D[i, j], 0.0, 0.0])  # fallback
                    kappa = dipole_orientation_factor(
                        dipoles[ki], dipoles[kj], R_vec)
                    J_ij *= abs(kappa)

            H[i, j] = H[j, i] = J_ij
    return H, keys


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python trp_extractor.py <PDB_ID> [--chain X] [--save]")
        sys.exit(1)

    pdb_id = sys.argv[1]
    chain_id = None
    save_csv = False
    if "--chain" in sys.argv:
        chain_id = sys.argv[sys.argv.index("--chain") + 1]
    if "--save" in sys.argv:
        save_csv = True

    pdb_text = fetch_pdb(pdb_id)
    if not pdb_text:
        sys.exit(1)

    centroids, dipoles = extract_trp_data(pdb_text, chain_id)
    if not centroids:
        print(f"[!] No Trp residues in {pdb_id}")
        sys.exit(0)

    print(f"\n{'='*60}")
    print(f"  TRYPTOPHAN NETWORK  -  {pdb_id.upper()}")
    print(f"  {len(centroids)} Trp residues")
    print(f"{'='*60}")
    for seq in sorted(centroids):
        c = centroids[seq]
        d = dipoles.get(seq, np.array([0,0,0]))
        print(f"  Trp-{seq:>4}:  ({c[0]:>8.3f}, {c[1]:>8.3f}, {c[2]:>8.3f})  mu=({d[0]:.3f},{d[1]:.3f},{d[2]:.3f})")

    D, keys = distance_matrix(centroids)
    coupled, relay = classify_pairs(D, keys)
    print(f"\n  Pairs: {len(coupled)} coupled, {len(relay)} relay")

    # Hamiltonian with dipole orientation
    H, h_keys = build_hamiltonian(D, keys, dipoles, centroids=centroids)
    print(f"\n  Hamiltonian [{len(h_keys)}x{len(h_keys)}] with dipole orientation:")
    for row in H:
        print(f"    [{', '.join(f'{x:>8.1f}' for x in row)}]")

    if save_csv:
        import csv
        out_dir = os.path.join(os.path.dirname(__file__), "..", "..", "data")
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f"{pdb_id.upper()}_trp_network.csv")
        with open(out_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Trp_i", "Trp_j", "Distance_A", "Regime", "Kappa"])
            for ki, kj, d in coupled + relay:
                regime = "coupled" if d < COUPLING_CUTOFF_ANGSTROM else "relay"
                kappa = 0.0
                if ki in dipoles and kj in dipoles:
                    Rv = centroids[kj] - centroids[ki]
                    kappa = dipole_orientation_factor(dipoles[ki], dipoles[kj], Rv)
                w.writerow([ki, kj, f"{d:.2f}", regime, f"{kappa:.3f}"])
        print(f"\n[+] Saved to {out_path}")
