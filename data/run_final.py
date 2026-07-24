import sys
sys.path.insert(0, 'src')
sys.path.insert(0, '.')
from pdb_tools.trp_extractor import fetch_pdb, extract_trp_data
from analysis.quantum_hamiltonian_engine import build_pdb_hamiltonian
import numpy as np

pdb_ids = ['1BL8','2BG9','6IRA','6AGF','4PE5','7KOX','6PM6','5C1M','3J5P','3GD8',
           '1SFC','3KG2','2A79','5GJV','5SVK','1U19','1JFF','3N2K','6HUO']
names = {'1BL8':'KcsA','2BG9':'nAChR','6IRA':'NMDA-A','6AGF':'Nav1.4','4PE5':'NMDA-B',
         '7KOX':'a7-nAChR','6PM6':'GlyR','5C1M':'mu-opioid','3J5P':'TRPV1','3GD8':'AQP4',
         '1SFC':'SNARE','3KG2':'AMPA','2A79':'Kv1.2','5GJV':'Cav1.1','5SVK':'P2X3',
         '1U19':'Rhod','1JFF':'Tubulin','3N2K':'Tub-mam','6HUO':'GABAA'}

print(f'{"PDB":>6s}  {"Name":>12s}  {"Trp":>3s}  {"S=2F":>8s}  {"Fcoh":>8s}  {"d_mean":>7s}')
print('-'*52)
for pid in pdb_ids:
    text = fetch_pdb(pid, cache=True)
    centroids, dipoles = extract_trp_data(text)
    keys = sorted(centroids.keys())
    n = len(keys)
    if n < 2:
        print(f'{pid}  {names.get(pid,"?"):>12s}  {n:>3d}  {"--":>8s}  {"--":>8s}  {"--":>7s}')
        continue
    H, _ = build_pdb_hamiltonian(centroids, dipoles)
    off = sum(abs(H[i,j].real) for i in range(n) for j in range(n) if i != j)
    diag = sum(abs(H[i,i].real) for i in range(n))
    F = off / (off + diag) if (off + diag) > 0 else 0
    S = 2.0 * F
    # mean distance
    all_d = []
    for i in range(n):
        for j in range(i+1, n):
            ri = centroids[keys[i]]
            rj = centroids[keys[j]]
            all_d.append(np.linalg.norm(np.array(ri)-np.array(rj)))
    d_mean = np.mean(all_d) if all_d else 0
    print(f'{pid:>6s}  {names.get(pid,"?"):>12s}  {n:>3d}  {S:>8.4f}  {F:>8.4f}  {d_mean:>7.1f}')
