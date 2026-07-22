import sys, numpy as np
sys.path.insert(0, 'src')
from pdb_tools.trp_extractor import fetch_pdb, extract_trp_data, distance_matrix, classify_pairs

pdb_ids = ['1BL8','6PV7','7TYO','6LQA','6CNO','7KOX','6J8J','1YAG','1JFF','3N2K',
           '6V2W','7SQS','6PLL','7S6M','6J4K','7EIX','6Z0C','7T6L','6N5B','7RJK']

names = {'1BL8':'KcsA K+','6PV7':'nAChR','7TYO':'NMDA-A','6LQA':'Nav1.4','6CNO':'NMDA-B',
         '7KOX':'a7-nAChR','6J8J':'Nav1.4b','1YAG':'SNARE','1JFF':'Tubulin','3N2K':'Tub-mam',
         '6V2W':'GABAA','7SQS':'AMPA','6PLL':'GlyR','7S6M':'Kv1.2','6J4K':'Cav1.1',
         '7EIX':'TRPV1','6Z0C':'P2X3','7T6L':'mu-opioid','6N5B':'Rhod','7RJK':'AQP4'}

print(f'{"PDB":>6s}  {"Protein":>10s}  {"Trp":>4s}  {"Coupled":>8s}  {"Relay":>6s}  {"MeanDist":>8s}')
print('-' * 55)
for pid in pdb_ids:
    text = fetch_pdb(pid, cache=True)
    if not text:
        print(f'{pid:>6s}  {"FAIL":>10s}')
        continue
    centroids, _ = extract_trp_data(text)
    if len(centroids) < 3:
        print(f'{pid:>6s}  {"<3 Trp":>10s}')
        continue
    keys = list(centroids.keys())
    D, _ = distance_matrix(centroids)
    coupled, relay = classify_pairs(D, keys)
    mean_d = np.mean(D[D > 0]) if np.sum(D > 0) else 0
    name = names.get(pid, '???')
    print(f'{pid:>6s}  {name:>10s}  {len(keys):>4d}  {len(coupled):>8d}  {len(relay):>6d}  {mean_d:>7.2f}A')
