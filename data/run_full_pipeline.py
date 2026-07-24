import sys, os, csv, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np

from pdb_tools.trp_extractor import fetch_pdb, extract_trp_coordinates, distance_matrix
from analysis.quantum_hamiltonian_engine import HamiltonianEngine

RESULTS_DIR = 'data'

# Corrected PDB list
pdb_ids = ['1BL8','2BG9','6IRA','6AGF','4PE5','7KOX','6PM6','5C1M','3J5P','3GD8',
           '1SFC','3KG2','2A79','5GJV','5SVK','1U19','1JFF','3N2K','6HUO']

names = {'1BL8':'KcsA','2BG9':'nAChR','6IRA':'NMDA-A','6AGF':'Nav1.4','4PE5':'NMDA-B',
         '7KOX':'a7-nAChR','6PM6':'GlyR','5C1M':'mu-opioid','3J5P':'TRPV1','3GD8':'AQP4',
         '1SFC':'SNARE','3KG2':'AMPA','2A79':'Kv1.2','5GJV':'Cav1.1','5SVK':'P2X3',
         '1U19':'Rhod','1JFF':'Tubulin','3N2K':'Tub-mam','6HUO':'GABAA'}

rows = []
for pid in pdb_ids:
    text = fetch_pdb(pid, cache=True)
    if not text:
        print(f'{pid} FAILED to fetch')
        continue
    centres = extract_trp_coordinates(text)
    n_trp = len(centres)
    if n_trp < 2:
        print(f'{pid} insufficient Trp ({n_trp})')
        rows.append({'pdb':pid,'name':names.get(pid,''),'n_trp':n_trp,'coupled':0,'relay':0,'mean_dist_A':0,'S_static':0,'error':'<2 Trp'})
        continue
    D, keys = distance_matrix(centres)
    n_pairs = n_trp * (n_trp - 1) // 2
    coupled = sum(1 for i in range(len(keys)) for j in range(i+1, len(keys)) if D[i,j] < 15.0)
    relay = n_pairs - coupled
    all_dists = [D[i,j] for i in range(len(keys)) for j in range(i+1, len(keys)) if D[i,j] > 0]
    mean_dist = float(np.mean(all_dists)) if all_dists else 0
    # Try HamiltonianEngine for KCKAS
    try:
        eng = HamiltonianEngine(pid)
        r = eng.scan_clock_drive(drives=[0, 0.25, 0.5, 0.75, 1.0])
        S_static = r['S_static']
        S_at_1 = r['S_at_clock_1']
    except Exception as e:
        S_static = 0
        S_at_1 = 0
        print(f'{pid} HEOM failed: {e}')
    rows.append({
        'pdb':pid,'name':names.get(pid,''),'n_trp':n_trp,'coupled':coupled,'relay':relay,
        'mean_dist_A':round(mean_dist,1),'S_static':round(S_static,4), 'S_at_1':round(S_at_1,4)
    })
    print(f'{pid:>6s} {names.get(pid,""):>10s} Trp={n_trp:2d} Coupled={coupled:2d} Relay={relay:3d} Dist={mean_dist:5.1f}A S0={S_static:.4f}')

# Save CSV
csv_path = os.path.join(RESULTS_DIR, 'corrected_pipeline_results.csv')
with open(csv_path, 'w', newline='') as f:
    if rows:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
print(f'\nSaved {len(rows)} results to {csv_path}')
