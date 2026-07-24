import sys
sys.path.insert(0, 'src')
sys.path.insert(0, '.')
from pdb_tools.trp_extractor import fetch_pdb, extract_trp_coordinates, distance_matrix
from analysis.quantum_hamiltonian_engine import HamiltonianEngine
import numpy as np

pdb_ids = ['1BL8','2BG9','6IRA','6AGF','4PE5','7KOX','6PM6','5C1M','3J5P','3GD8',
           '1SFC','3KG2','2A79','5GJV','5SVK','1U19','1JFF','3N2K','6HUO']
names = {'1BL8':'KcsA','2BG9':'nAChR','6IRA':'NMDA-A','6AGF':'Nav1.4','4PE5':'NMDA-B',
         '7KOX':'a7-nAChR','6PM6':'GlyR','5C1M':'mu-opioid','3J5P':'TRPV1','3GD8':'AQP4',
         '1SFC':'SNARE','3KG2':'AMPA','2A79':'Kv1.2','5GJV':'Cav1.1','5SVK':'P2X3',
         '1U19':'Rhod','1JFF':'Tubulin','3N2K':'Tub-mam','6HUO':'GABAA'}

results = []
for pid in pdb_ids:
    text = fetch_pdb(pid, cache=True)
    if not text:
        results.append(f'{pid} FAILED')
        continue
    centres = extract_trp_coordinates(text)
    n_trp = len(centres)
    if n_trp < 2:
        results.append(f'{pid}  {names[pid]:>10s}  {n_trp:>2d}  {"--":>8s}  {"--":>8s}  {"--":>8s}')
        continue
    D, keys = distance_matrix(centres)
    n_pairs = n_trp * (n_trp - 1) // 2
    coupled = sum(1 for i in range(n_trp) for j in range(i+1, n_trp) if D[i,j] < 15.0)
    relay = n_pairs - coupled
    all_d = [D[i,j] for i in range(n_trp) for j in range(i+1, n_trp) if D[i,j] > 0]
    mean_d = float(np.mean(all_d)) if all_d else 0
    S_static = 0.0
    try:
        eng = HamiltonianEngine(pid)
        r = eng.scan_clock_drive(drives=[0, 0.25])
        S_static = r['S_static']
    except:
        pass
    results.append(f'{pid}  {names[pid]:>10s}  {n_trp:>2d}  {coupled:>2d}  {relay:>3d}  {mean_d:6.1f}A  S0={S_static:.4f}')

print(f'\n{"PDB":>6s}  {"Protein":>10s}  {"Trp":>2s}  {"Cpl":>2s}  {"Rel":>3s}  {"Dist":>8s}  {"S0":>8s}')
print('-'*55)
for r in results:
    print(r)
