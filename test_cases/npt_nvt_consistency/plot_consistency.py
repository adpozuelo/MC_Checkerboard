#!/usr/bin/env python3
"""
plot_consistency.py: Generates publication-quality comparison plots
of density and energy distributions between NpT and NVT ensembles.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_virial_p(fpath, is_hs=False):
    steps, p_vir = [], []
    col = 1 if is_hs else 2
    if os.path.exists(fpath):
        with open(fpath) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= (col + 1) and parts[0].isdigit():
                    try:
                        steps.append(float(parts[0]))
                        p_vir.append(float(parts[col]))
                    except ValueError:
                        pass
    return np.array(steps), np.array(p_vir)

def plot_system(system_dir, system_name):
    npt_file = os.path.join(system_dir, 'npt', 'run-data.dat')
    nvt_file = os.path.join(system_dir, 'nvt', 'run-data.dat')
    
    if not os.path.exists(npt_file) or not os.path.exists(nvt_file):
        return False
        
    try:
        d_npt = np.loadtxt(npt_file)
        d_nvt = np.loadtxt(nvt_file)
        if d_npt.ndim < 2 or d_nvt.ndim < 2 or len(d_npt) == 0 or len(d_nvt) == 0:
            return False
    except Exception:
        return False
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # 1. Density evolution
    ax1.plot(d_npt[:, 0], d_npt[:, 2], label=r'NpT Ensemble $\rho^*(t)$', color='#1f77b4', lw=1.5)
    rho_nvt_val = np.mean(d_nvt[:, 2])
    ax1.axhline(rho_nvt_val, color='#d62728', linestyle='--', lw=2, label=rf'NVT Fixed $\langle \rho^* \rangle = {rho_nvt_val:.4f}$')
    ax1.set_xlabel('Monte Carlo Sweeps', fontsize=12)
    ax1.set_ylabel(r'Reduced Density $\rho^*$', fontsize=12)
    ax1.set_title(f'{system_name.upper()}: Density Evolution', fontsize=13, fontweight='bold')
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(frameon=True, facecolor='white')
    
    # 2. Virial Pressure evolution
    is_hs = (system_name.lower() == 'hs')
    s_npt, p_npt = parse_virial_p(os.path.join(system_dir, 'npt', 'output.out'), is_hs=is_hs)
    s_nvt, p_nvt = parse_virial_p(os.path.join(system_dir, 'nvt', 'output.out'), is_hs=is_hs)
    
    if len(s_npt) > 0:
        ax2.plot(s_npt, p_npt, label=r'NpT $P_{\rm vir}^*$', color='#2ca02c', alpha=0.7, lw=1.2)
    if len(s_nvt) > 0:
        ax2.plot(s_nvt, p_nvt, label=r'NVT $P_{\rm vir}^*$', color='#ff7f0e', alpha=0.7, lw=1.2)
        
    if system_name.lower() == 'hs':
        eta_nvt = (np.pi / 6.0) * rho_nvt_val
        p_ref = rho_nvt_val * (1.0 + eta_nvt + eta_nvt**2 - eta_nvt**3) / ((1.0 - eta_nvt)**3)
        ref_label = f'Carnahan-Starling ({p_ref:.4f})'
    elif system_name.lower() == 'lj':
        p_ref = 0.60
        ref_label = f'Target $P^* = {p_ref:.4f}$'
    elif system_name.lower() == 'ssp':
        p_ref = 0.05
        ref_label = f'Target $P^* = {p_ref:.4f}$'
    else:
        p_ref = None

    if p_ref is not None:
        ax2.axhline(p_ref, color='#d62728', linestyle='--', lw=2, label=ref_label)

    ax2.set_xlabel('Monte Carlo Sweeps', fontsize=12)
    ax2.set_ylabel(r'Virial Pressure $P^* = P\sigma^3/k_B T$', fontsize=12)
    ax2.set_title(f'{system_name.upper()}: Virial Pressure Consistency', fontsize=13, fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(frameon=True, facecolor='white')
    
    plt.tight_layout()
    out_png = os.path.join(system_dir, f'consistency_{system_name.lower()}.png')
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"Plot saved to: {out_png}")
    return True

def main():
    base_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    for s in ['hs', 'lj', 'ssp']:
        s_dir = os.path.join(base_dir, s)
        if os.path.isdir(s_dir):
            plot_system(s_dir, s)

if __name__ == '__main__':
    main()
