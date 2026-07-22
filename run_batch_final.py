import sys, gc, csv
sys.path.insert(0, 'src')
from src.v2.multiscale_pipeline import run_pipeline

targets = ['1BL8','6PV7','7TYO','6LQA','6CNO','7KOX','6J8J','1YAG','1JFF','3N2K',
           '6V2W','7SQS','6PLL','7S6M','7EIX','7T6L','6N5B','7RJK']

results = {}
for pid in targets:
    print(f'\n=== {pid} ===')
    r = run_pipeline(pid, verbose=False, skip_heom=True)
    results[pid] = r
    s = r.get('kckas_score', 0)
    ac = r.get('A_crit', 0)
    ps = r.get('p_success', 0)
    print(f'  KCKAS={s:.4f}  A_crit={ac:.2f}  P_succ={ps:.4%}')
    gc.collect()

keys = ['pdb_id','n_trp','J_max_cm','J_mean_cm','tau_coherence_fs','f_coh','kckas_score',
        'A_crit','S_clock_max','clock_rescues','p_hit','p_success','capacity_bits','e_bit_J','landauer_ratio']
with open('data/pipeline_results.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    for pid in targets:
        w.writerow({k: results[pid].get(k, '') for k in keys})
print('\nAll complete, CSV saved')
