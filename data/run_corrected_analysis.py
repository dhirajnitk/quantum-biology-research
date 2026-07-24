import sys, numpy as np
sys.path.insert(0, 'src')
from pdb_tools.trp_extractor import fetch_pdb, extract_trp_data, distance_matrix, classify_pairs

pdb_ids = ['1BL8','2BG9','6IRA','6AGF','4PE5','7KOX','6PM6','5C1M','3J5P','3GD8',
           '1SFC','3KG2','2A79','5GJV','5SVK','1U19','1JFF','3N2K','6HUO']

names = {'1BL8':'KcsA K+','2BG9':'nAChR','6IRA':'NMDA-A','6AGF':'Nav1.4','4PE5':'NMDA-B',
         '7KOX':'a7-nAChR','6PM6':'GlyR','5C1M':'mu-opioid','3J5P':'TRPV1','3GD8':'AQP4',
         '1SFC':'SNARE','3KG2':'AMPA','2A79':'Kv1.2','5GJV':'Cav1.1','5SVK':'P2X3',
         '1U19':'Rhod','1JFF':'Tubulin','3N2K':'Tub-mam','6HUO':'GABAA'}

print(f'{"PDB":>6s}  {"Protein":>10s}  {"Trp":>4s}  {"Coupled":>8s}  {"Relay":>6s}  {"MeanDist":>8s}')
print('-' * 55)
for pid in pdb_ids:
    text = fetch_pdb(pid, cache=True)
    if not text:
        print(f'{pid:>6s}  {"FAIL":>10s}')
        continue
    centroids, _ = extract_trp_data(text)
    if len(centroids) < 2:
        print(f'{pid:>6s}  {"SKIP":>10s}  {"<2 Trp":>4s}')
        continue
    keys = list(centroids.keys())
    D, _ = distance_matrix(centroids)
    coupled, relay = classify_pairs(D, keys)
    mean_d = np.mean(D[D > 0]) if np.sum(D > 0) else 0
    name = names.get(pid, '???')
    print(f'{pid:>6s}  {name:>10s}  {len(keys):>4d}  {len(coupled):>8d}  {len(relay):>6d}  {mean_d:>7.2f}A')
