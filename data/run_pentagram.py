import sys
sys.path.insert(0, 'src')
sys.path.insert(0, '.')
sys.path.insert(0, 'src/analysis')
import numpy as np
from pdb_contextuality import select_trp_pentagram
from pdb_tools.trp_extractor import fetch_pdb, extract_trp_data
from analysis.quantum_hamiltonian_engine import build_pdb_hamiltonian

pdb_ids = ['1BL8','2BG9','6IRA','6AGF','4PE5','7KOX','6PM6','5C1M','3J5P','3GD8',
           '1SFC','3KG2','2A79','5GJV','5SVK','1U19','1JFF','3N2K','6HUO']
names = {'1BL8':'KcsA','2BG9':'nAChR','6IRA':'NMDA-A','6AGF':'Nav1.4','4PE5':'NMDA-B',
         '7KOX':'a7-nAChR','6PM6':'GlyR','5C1M':'mu-opioid','3J5P':'TRPV1','3GD8':'AQP4',
         '1SFC':'SNARE','3KG2':'AMPA','2A79':'Kv1.2','5GJV':'Cav1.1','5SVK':'P2X3',
         '1U19':'Rhod','1JFF':'Tubulin','3N2K':'Tub-mam','6HUO':'GABAA'}

print(f'{"PDB":>6s}  {"Name":>12s}  {"Trp":>3s}  {"S=2F":>8s}  {"Fcoh":>8s}  {"d_mean(A)":>9s}')
print('-'*52)
for pid in pdb_ids:
    text = fetch_pdb(pid, cache=True)
    centroids, dipoles = extract_trp_data(text)
    # Use fixed selection algorithm (tightest 5-residue cluster)
    result = select_trp_pentagram(centroids, n_nodes=5)
    if result is None:
        print(f'{pid}  {names.get(pid,"?"):>12s}  {"?":>3s}  {"--":>8s}  {"--":>8s}  {"--":>9s}')
        continue
    sel_coords, sel_keys = result
    # Subset centroids to the 5 selected residues
    selected_centroids = {k: centroids[k] for k in sel_keys}
    selected_dipoles = {k: dipoles[k] for k in sel_keys}
    n_sel = len(sel_keys)
    # Build Hamiltonian and compute F_coh
    H, _ = build_pdb_hamiltonian(selected_centroids, selected_dipoles)
    off = sum(abs(H[i,j].real) for i in range(n_sel) for j in range(n_sel) if i != j)
    diag = sum(abs(H[i,i].real) for i in range(n_sel))
    F = off / (off + diag) if (off + diag) > 0 else 0
    S = 2.0 * F
    # Mean distance among selected residues
    all_d = []
    for i in range(n_sel):
        for j in range(i+1, n_sel):
            ri, rj = selected_centroids[sel_keys[i]], selected_centroids[sel_keys[j]]
            all_d.append(np.linalg.norm(np.array(ri)-np.array(rj)))
    d_mean = np.mean(all_d) if all_d else 0
    print(f'{pid:>6s}  {names.get(pid,"?"):>12s}  {n_sel:>3d}  {S:>8.4f}  {F:>8.4f}  {d_mean:>9.1f}')
