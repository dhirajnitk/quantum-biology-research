"""
ray_tracing_3d.py — 3D Monte Carlo optical ray-tracing for CCO path loss.

Simulates biophoton propagation from Trp sources in the presynaptic active
zone membrane to CCO targets in the inner mitochondrial membrane, accounting
for solid-angle attenuation, membrane waveguiding, and the intermembrane gap.
"""

import numpy as np
from numpy.random import uniform, exponential


def simulate_photon_burst(
    n_photons=200000,
    az_radius_nm=200.0,       # active zone radius
    cleft_nm=20.0,             # synaptic cleft width
    ims_gap_nm=15.0,           # intermembrane space gap
    membrane_n=1.45,           # lipid refractive index
    cyto_n=1.33,               # cytoplasm refractive index
    cco_density=50,            # CCO complexes per um^2
    cco_sigma_m2=8.0e-21,      # CCO absorption cross-section
    membrane_loss=0.002,       # loss coeff in membrane (nm^-1)
    seed=42,
):
    """Run 3D ray-tracing from active zone to mitochondrial CCO.

    Geometry:
        z=0: presynaptic membrane (Trp sources)
        z=cleft_nm: outer mitochondrial membrane
        z=cleft_nm + ims_gap_nm: inner mitochondrial membrane (CCO)

    Returns
    -------
    dict with p_hit, hits, n_photons, details
    """
    np.random.seed(seed)

    # CCO targets on inner mitochondrial membrane
    n_cco = int(cco_density * np.pi * az_radius_nm**2 * 1e-6)  # per active zone
    cco_x = uniform(-az_radius_nm, az_radius_nm, n_cco)
    cco_y = uniform(-az_radius_nm, az_radius_nm, n_cco)
    # Keep only those within the active zone disk
    mask = cco_x**2 + cco_y**2 < az_radius_nm**2
    cco_x = cco_x[mask]
    cco_y = cco_y[mask]
    n_cco = len(cco_x)
    r_cco = np.sqrt(cco_sigma_m2 * 1e18 / np.pi)  # effective radius in nm

    target_z = cleft_nm + ims_gap_nm

    # Trp source positions (uniform in active zone disk)
    theta = uniform(0, 2*np.pi, n_photons)
    r = az_radius_nm * np.sqrt(uniform(0, 1, n_photons))
    x0 = r * np.cos(theta)
    y0 = r * np.sin(theta)

    hits = 0
    for i in range(n_photons):
        # Random direction (isotropic into forward hemisphere)
        phi = uniform(0, 2*np.pi)
        cos_theta = uniform(0, 1)
        dx = np.sin(np.arccos(cos_theta)) * np.cos(phi)
        dy = np.sin(np.arccos(cos_theta)) * np.sin(phi)
        dz = cos_theta

        # Step through membrane and IMS
        x, y, z = x0[i], y0[i], 0.0
        absorbed = False

        while z < target_z and z >= 0:
            step = 1.0  # nm

            # In lipid membrane (z < cleft_nm): lower loss, waveguide
            # In IMS (z >= cleft_nm): higher loss, aqueous
            if z < cleft_nm:
                loss = membrane_loss
            else:
                loss = 0.01  # higher loss in aqueous IMS

            x += dx * step
            y += dy * step
            z += dz * step

            # Absorption by medium
            if uniform() < loss * step:
                absorbed = True
                break

            # Scattering
            scatter_prob = 0.001 if z < cleft_nm else 0.01
            if uniform() < scatter_prob * step:
                phi = uniform(0, 2*np.pi)
                ct = uniform(-1, 1)
                st = np.sqrt(1 - ct**2)
                dx = st * np.cos(phi)
                dy = st * np.sin(phi)
                dz = ct if z < cleft_nm else abs(ct)  # bias forward in IMS

        if absorbed or z < target_z:
            continue

        # Check CCO absorption
        for j in range(n_cco):
            if (x - cco_x[j])**2 + (y - cco_y[j])**2 < r_cco**2:
                hits += 1
                break

    p_hit = hits / n_photons if n_photons > 0 else 0.0
    return {
        'p_hit': p_hit,
        'hits': hits,
        'n_photons': n_photons,
        'n_cco': n_cco,
        'r_cco': r_cco,
    }


def ensemble_from_raytrace(p_hit, N=5000):
    """Compute ensemble gating reliability from ray-traced p_hit."""
    return 1 - (1 - p_hit) ** N


if __name__ == '__main__':
    print("=== 3D Ray-Tracing: CCO Path Loss ===\n")

    # Default scenario
    result = simulate_photon_burst(n_photons=50000)
    p = result['p_hit']
    print(f"Photons: {result['n_photons']}")
    print(f"CCO targets: {result['n_cco']}")
    print(f"CCO effective radius: {result['r_cco']:.3f} nm")
    print(f"Hits: {result['hits']}")
    print(f"p_hit (ray-trace) = {p:.6e}")
    print(f"Paper p = 7.84e-4")
    print(f"Ensemble P_success(N=5000): ray-trace = {ensemble_from_raytrace(p):.4%}, "
          f"paper = {ensemble_from_raytrace(7.84e-4):.4%}")

    # Sweep IMS gap
    print("\n--- IMS Gap Sweep ---")
    for gap in [5, 10, 15, 20, 30]:
        r = simulate_photon_burst(n_photons=20000, ims_gap_nm=gap)
        print(f"  IMS gap={gap}nm: p={r['p_hit']:.2e}, P_success(5000)={ensemble_from_raytrace(r['p_hit']):.3%}")
