import numpy as np

def kinetic_bernoulli_mc(N=5000, p=5.9e-4, trials=200000, threshold=1):
    successes = np.random.binomial(n=N, p=p, size=trials)
    gating = (successes >= threshold).astype(float)
    return np.mean(gating), np.std(gating) / np.sqrt(trials)

def kinetic_bernoulli_robustness(trials=200000):
    p_fixed = 5.9e-4
    N_fixed = 5000
    results = {}

    s = np.random.binomial(n=N_fixed, p=p_fixed, size=trials)
    results['ideal'] = np.mean(s > 0)

    N_pois = np.random.poisson(N_fixed, size=trials)
    alpha_p, beta_p = 590, 999410
    p_beta = np.random.beta(alpha_p, beta_p, size=trials)
    s = np.random.binomial(n=N_pois, p=p_beta)
    results['biological_noise'] = np.mean(s > 0)

    s = np.random.binomial(n=N_pois, p=p_beta)
    results['threshold_k2'] = np.mean(s >= 2)

    N_dmg = np.random.poisson(4000, size=trials)
    s = np.random.binomial(n=N_dmg, p=p_fixed)
    results['protein_damage'] = np.mean(s > 0)

    p_stress = np.random.beta(531, 899469, size=trials)
    s = np.random.binomial(n=N_fixed, p=p_stress)
    results['oxidative_stress'] = np.mean(s > 0)

    return results

def required_ensemble_for_target(targets=None, p=5.9e-4, trials=200000):
    if targets is None:
        targets = [0.50, 0.90, 0.95, 0.99, 0.999]
    Ns = np.logspace(1, 5, 200).astype(int)
    results = {}
    for target in targets:
        ideal = next(N for N in Ns if 1 - (1-p)**N >= target)
        poisson = next(N for N in Ns if np.mean(
            np.random.binomial(n=np.random.poisson(N, size=trials), p=p) > 0) >= target)
        damage = next(N for N in Ns if np.mean(
            np.random.binomial(n=np.random.poisson(int(N*0.8), size=trials), p=p) > 0) >= target)
        results[f'{target*100:.1f}%'] = {'ideal': ideal, 'poisson_N': poisson, 'damage_20pct': damage}
    return results

def markov_cco_ensemble(N=5000, p_photo=5.9e-4, trials=10000, t_max=100.0, dt=0.1):
    k_relax = 0.1
    k_pump = 1.0
    n_steps = int(t_max / dt)
    gates = 0
    for _ in range(trials):
        n_photons = np.random.binomial(N, p_photo)
        state = 0
        gated = False
        for _ in range(n_steps):
            if state == 0 and n_photons > 0:
                state = 1
                n_photons = 0
            elif state == 1:
                r = np.random.uniform()
                if r < k_relax * dt:
                    state = 0
                elif r < (k_relax + k_pump) * dt:
                    gated = True
                    break
        if gated:
            gates += 1
    return gates / trials

def markov_burst_gamma(N=5000, n_bursts=10, trials=5000, gap_ms=25):
    dt_ps = 0.1
    gap_steps = int(gap_ms * 1000 / dt_ps)
    gates = 0
    for _ in range(trials):
        gated = False
        for _ in range(n_bursts):
            if np.random.binomial(N, 5.9e-4) > 0:
                gated = True
                break
        if gated:
            gates += 1
    return gates / trials

if __name__ == '__main__':
    print('=== Kinetic Bernoulli MC ===')
    r, e = kinetic_bernoulli_mc()
    print(f'Ideal (N=5000, p=5.9e-4):  {r*100:.2f}%  +/- {e*100:.2f}%')

    print('\n=== Robustness Analysis ===')
    rob = kinetic_bernoulli_robustness()
    for k, v in rob.items():
        print(f'  {k}: {v*100:.2f}%')

    print('\n=== Required Ensemble Sizes ===')
    req = required_ensemble_for_target()
    for target, vals in req.items():
        print(f'  {target}: ideal={vals["ideal"]}, Poisson_N={vals["poisson_N"]}, damage={vals["damage_20pct"]}')

    print('\n=== Markov Chain CCO ===')
    for N in [1000, 5000, 10000]:
        g = markov_cco_ensemble(N=N, trials=5000)
        print(f'  CCO gating N={N}:  {g*100:.2f}%')

    print('\n  Gamma burst gating (40 Hz):')
    for N in [1000, 5000, 10000]:
        g = markov_burst_gamma(N=N)
        print(f'  N={N}:  {g*100:.2f}%')

    print('\n=== Dielectric Perturbation MC ===')
    trials = 200000
    N = 5000
    p0 = 5.9e-4
    P_fixed = np.mean(np.random.binomial(N, p0, size=trials) > 0)
    eps = np.random.uniform(1.8, 2.5, size=trials)
    p = p0 * np.sqrt(2.0 / eps)
    P_pert = np.mean(1 - (1 - p) ** N)
    eps2 = np.random.uniform(1.5, 3.0, size=(trials, N))
    p2 = p0 * np.sqrt(2.0 / eps2)
    P_extreme = np.mean(1 - np.prod(1 - p2, axis=1))
    print(f'  Fixed epsilon=2.0:             {P_fixed*100:.2f}%')
    print(f'  Perturbed [1.8,2.5] (uniform): {P_pert*100:.2f}%')
    print(f'  Per-core [1.5,3.0] (worst):    {P_extreme*100:.2f}%')
