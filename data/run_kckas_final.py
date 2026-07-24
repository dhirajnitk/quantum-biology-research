import sys, os
sys.path.insert(0, 'src')
sys.path.insert(0, '.')
from pdb_tools.trp_extractor import fetch_pdb, extract_trp_data
from analysis.quantum_hamiltonian_engine import build_pdb_hamiltonian, compute_kckas_sdp
import numpy as np

pdb_ids = ['1BL8','2BG9','6IRA','6AGF','4PE5','7KOX','6PM6','5C1M','3J5P','3GD8',
           '1SFC','3KG2','2A79','5GJV','5SVK','1U19','1JFF','3N2K','6HUO']
names = {'1BL8':'KcsA','2BG9':'nAChR','6IRA':'NMDA-A','6AGF':'Nav1.4','4PE5':'NMDA-B',
         '7KOX':'a7-nAChR','6PM6':'GlyR','5C1M':'mu-opioid','3J5P':'TRPV1','3GD8':'AQP4',
         '1SFC':'SNARE','3KG2':'AMPA','2A79':'Kv1.2','5GJV':'Cav1.1','5SVK':'P2X3',
         '1U19':'Rhod','1JFF':'Tubulin','3N2K':'Tub-mam','6HUO':'GABAA'}

print(f'{"PDB":>6s}  {"Name":>10s}  {"Trp":>3s}  {"S":>8s}  {"Fcoh":>8s}')
print('-'*42)
for pid in pdb_ids:
    text = fetch_pdb(pid, cache=True)
    if not text:
        print(f'{pid}  FAIL')
        continue
    centroids, dipoles = extract_trp_data(text)
    keys = sorted(centroids.keys())
    n = len(keys)
    if n < 2:
        print(f'{pid}  {names.get(pid,"?"):>10s}  {n:>3d}  {"--":>8s}  {"--":>8s}')
        continue
    # Build Hamiltonian with site energies
    H, _ = build_pdb_hamiltonian(centroids, dipoles)
    S = compute_kckas_sdp(H)
    # F_coh from actual matrix
    off = sum(abs(H[i,j].real) for i in range(n) for j in range(n) if i != j)
    diag = sum(abs(H[i,i].real) for i in range(n))
    F = off / (off + diag) if (off + diag) > 0 else 0
    print(f'{pid:>6s}  {names.get(pid,"?"):>10s}  {n:>3d}  {S:>8.4f}  {F:>8.4f}')
