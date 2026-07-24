import sys, os
sys.path.insert(0, 'src')
sys.path.insert(0, '.')
from analysis.quantum_hamiltonian_engine import HamiltonianEngine
import numpy as np

pdb_ids = ['1BL8','2BG9','6IRA','6AGF','4PE5','7KOX','6PM6','5C1M','3J5P','3GD8',
           '1SFC','3KG2','2A79','5GJV','5SVK','1U19','1JFF','3N2K','6HUO']
names = {'1BL8':'KcsA','2BG9':'nAChR','6IRA':'NMDA-A','6AGF':'Nav1.4','4PE5':'NMDA-B',
         '7KOX':'a7-nAChR','6PM6':'GlyR','5C1M':'mu-opioid','3J5P':'TRPV1','3GD8':'AQP4',
         '1SFC':'SNARE','3KG2':'AMPA','2A79':'Kv1.2','5GJV':'Cav1.1','5SVK':'P2X3',
         '1U19':'Rhod','1JFF':'Tubulin','3N2K':'Tub-mam','6HUO':'GABAA'}

print(f'{"PDB":>6s}  {"Name":>12s}  {"Trp":>3s}  {"S_static":>10s}  {"Status":>12s}')
print('-'*48)
for pid in pdb_ids:
    try:
        eng = HamiltonianEngine(pid)
        r = eng.scan_clock_drive(drives=[0])
        S = r['S_static']
        print(f'{pid:>6s}  {names.get(pid,"?"):>12s}  {eng.n_trp:>3d}  {S:>10.4f}  {"OK":>12s}')
    except Exception as e:
        msg = str(e)[:40] if str(e) else 'unknown'
        print(f'{pid:>6s}  {names.get(pid,"?"):>12s}  {"?":>3s}  {"FAIL":>10s}  {msg:>12s}')
