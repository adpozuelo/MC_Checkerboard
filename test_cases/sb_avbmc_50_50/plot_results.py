#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_run_data(fname):
    steps = []
    energies = []
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 8:
                try:
                    s = int(parts[0])
                    e = float(parts[7])
                    steps.append(s)
                    energies.append(e)
                except ValueError:
                    continue
    return np.array(steps), np.array(energies)

def main():
    steps, energy = parse_run_data('run-data.dat')
    print('=== Simulation Data Analysis ===')
    print(f'Total steps executed: {steps[-1]}')
    print(f'Initial energy per particle: {energy[0]:.4f} kBT')
    print(f'Final energy per particle: {energy[-1]:.4f} kBT')

    # Equilibrium window: steps >= 30000
    mask_plateau = steps >= 30000
    e_plat = energy[mask_plateau]
    print(f'Plateau (steps >= 30k) mean energy per particle: {np.mean(e_plat):.4f} +/- {np.std(e_plat):.4f} kBT')

    # Clusevol data
    cl_steps, n_cl, max_sz, pct_cl = [], [], [], []
    with open('clusevol_mc.dat', 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    cl_steps.append(int(parts[0]))
                    n_cl.append(int(parts[1]))
                    max_sz.append(int(parts[2]))
                    pct_cl.append(float(parts[4]))
                except ValueError:
                    continue
    cl_steps = np.array(cl_steps)
    n_cl = np.array(n_cl)
    max_sz = np.array(max_sz)
    pct_cl = np.array(pct_cl)

    mask_cl = cl_steps >= 30000
    print(f'Plateau mean largest cluster size: {np.mean(max_sz[mask_cl]):.1f} +/- {np.std(max_sz[mask_cl]):.1f} particles')
    print(f'Plateau mean % clustered: {np.mean(pct_cl[mask_cl]):.2f}% +/- {np.std(pct_cl[mask_cl]):.2f}%')
    print(f'Plateau mean number of clusters: {np.mean(n_cl[mask_cl]):.1f} +/- {np.std(n_cl[mask_cl]):.1f}')

    # Create publication quality plot
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    # 1. Energy evolution
    ax1.plot(steps, energy, color='#1f77b4', lw=1.2, label=r'$E_{\rm tot}/N$')
    ax1.axhline(np.mean(e_plat), color='red', linestyle='--', lw=1.0, label=f'Plateau: {np.mean(e_plat):.2f} $k_B T$')
    ax1.set_ylabel('Energy per particle ($k_B T$)', fontsize=11)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc='lower right', frameon=True)
    ax1.set_title('Sanchez-Burgos 50:50 Mixture (AVBMC + Tables)', fontsize=13, fontweight='bold')

    # 2. Largest cluster size & % clustered
    ax2.plot(cl_steps, max_sz, color='#2ca02c', lw=1.5, label='Max Cluster Size')
    ax2.axhline(2000, color='gray', linestyle=':', label='Total System (N=2000)')
    ax2.set_ylabel('Max Cluster Size (particles)', fontsize=11, color='#2ca02c')
    ax2.tick_params(axis='y', labelcolor='#2ca02c')
    ax2.grid(True, linestyle=':', alpha=0.6)

    ax2_r = ax2.twinx()
    ax2_r.plot(cl_steps, pct_cl, color='#ff7f0e', lw=1.2, linestyle='--', label='% Clustered')
    ax2_r.set_ylabel('Clustered Fraction (%)', fontsize=11, color='#ff7f0e')
    ax2_r.tick_params(axis='y', labelcolor='#ff7f0e')
    ax2_r.set_ylim(0, 105)

    # 3. Cluster Count
    ax3.plot(cl_steps, n_cl, color='#9467bd', lw=1.2)
    ax3.set_ylabel('Cluster Count', fontsize=11)
    ax3.set_xlabel('Monte Carlo Sweeps', fontsize=11)
    ax3.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig('sb_avbmc_evolution.png', dpi=300)
    print("Saved plot to sb_avbmc_evolution.png")

if __name__ == '__main__':
    main()
