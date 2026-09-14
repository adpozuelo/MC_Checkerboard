#!/usr/bin/env python3
"""
plot_consistency_grouped.py
Generates publication-quality thermodynamic consistency figures in PDF and PNG.
Groups figures by system:
- Hard Spheres (HS)
- Lennard-Jones (LJ)
- Site-Site Patchy (SSP)
For each system, density (top) and energy/virial pressure (bottom) are stacked
vertically sharing the x-axis for maximum readability.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Publication styling
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9.5,
    'figure.titlesize': 14,
    'lines.linewidth': 1.8,
    'figure.dpi': 300
})

def parse_hs_p(fpath):
    steps, p_vir = [], []
    if os.path.exists(fpath):
        with open(fpath) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 8 and parts[0].isdigit():
                    try:
                        steps.append(float(parts[0]))
                        p_vir.append(float(parts[1]))
                    except ValueError:
                        pass
    return np.array(steps), np.array(p_vir)

def load_data(base_dir):
    data = {}
    for sys_name in ['hs', 'lj', 'ssp']:
        sys_dir = os.path.join(base_dir, sys_name)
        npt_file = os.path.join(sys_dir, 'npt', 'run-data.dat')
        nvt_file = os.path.join(sys_dir, 'nvt', 'run-data.dat')
        d_npt = np.loadtxt(npt_file)
        d_nvt = np.loadtxt(nvt_file)
        
        sys_dict = {
            'd_npt': d_npt,
            'd_nvt': d_nvt,
        }
        
        if sys_name == 'hs':
            s_npt, p_npt = parse_hs_p(os.path.join(sys_dir, 'npt', 'output.out'))
            s_nvt, p_nvt = parse_hs_p(os.path.join(sys_dir, 'nvt', 'output.out'))
            sys_dict['s_npt'] = s_npt
            sys_dict['p_npt'] = p_npt
            sys_dict['s_nvt'] = s_nvt
            sys_dict['p_nvt'] = p_nvt
        data[sys_name] = sys_dict
    return data

def generate_grouped_figure(data, out_dir):
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.2), sharex='col')
    plt.subplots_adjust(hspace=0.08, wspace=0.28, top=0.92, bottom=0.10, left=0.07, right=0.98)
    
    # ------------------ COLUMN 1: HARD SPHERES (HS) ------------------
    d_hs = data['hs']
    ax_top = axes[0, 0]
    ax_bot = axes[1, 0]
    
    # Top: Density
    ax_top.plot(d_hs['d_npt'][:, 0], d_hs['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
    ax_top.axhline(np.mean(d_hs['d_nvt'][:, 2]), color='#d62728', linestyle='--', lw=2.0,
                   label=r'$NVT$ $\langle \rho^* \rangle = 0.3364$')
    ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=12, fontweight='bold')
    ax_top.set_title('(a) Hard Spheres (HS)', fontsize=13, fontweight='bold', pad=8)
    ax_top.grid(True, linestyle=':', alpha=0.6)
    ax_top.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax_top.set_ylim(0.295, 0.345)
    
    # Bottom: Virial Pressure
    s_npt, p_npt = d_hs['s_npt'], d_hs['p_npt']
    s_nvt, p_nvt = d_hs['s_nvt'], d_hs['p_nvt']
    eta_nvt = (np.pi / 6.0) * np.mean(d_hs['d_nvt'][:, 2])
    p_cs = np.mean(d_hs['d_nvt'][:, 2]) * (1.0 + eta_nvt + eta_nvt**2 - eta_nvt**3) / ((1.0 - eta_nvt)**3)
    
    if len(s_npt) > 0:
        ax_bot.plot(s_npt, p_npt, label=r'$NpT$ $P_{\rm vir}^*$', color='#2ca02c', alpha=0.75, lw=1.3)
    if len(s_nvt) > 0:
        ax_bot.plot(s_nvt, p_nvt, label=r'$NVT$ $P_{\rm vir}^*$', color='#ff7f0e', alpha=0.85, lw=1.3)
    ax_bot.axhline(p_cs, color='#d62728', linestyle='--', lw=2.0, label=f'CS EOS ({p_cs:.4f})')
    ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=12, fontweight='bold')
    ax_bot.set_ylabel(r'Virial Pressure $P\sigma^3/k_B T$', fontsize=12, fontweight='bold')
    ax_bot.grid(True, linestyle=':', alpha=0.6)
    ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax_bot.set_ylim(0.55, 0.80)
    
    # ------------------ COLUMN 2: LENNARD-JONES (LJ) ------------------
    d_lj = data['lj']
    ax_top = axes[0, 1]
    ax_bot = axes[1, 1]
    
    # Top: Density
    ax_top.plot(d_lj['d_npt'][:, 0], d_lj['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
    ax_top.axhline(np.mean(d_lj['d_nvt'][:, 2]), color='#d62728', linestyle='--', lw=2.0,
                   label=r'$NVT$ $\langle \rho^* \rangle = 0.4703$')
    ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=12, fontweight='bold')
    ax_top.set_title('(b) Lennard-Jones (LJ)', fontsize=13, fontweight='bold', pad=8)
    ax_top.grid(True, linestyle=':', alpha=0.6)
    ax_top.legend(loc='upper right', frameon=True, framealpha=0.9)
    ax_top.set_ylim(0.445, 0.585)
    
    # Bottom: Energy
    ax_bot.plot(d_lj['d_npt'][:, 0], d_lj['d_npt'][:, 7], label=r'$NpT$ $\langle U/N \rangle$', color='#2ca02c', alpha=0.85, lw=1.3)
    ax_bot.plot(d_lj['d_nvt'][:, 0], d_lj['d_nvt'][:, 7], label=r'$NVT$ $\langle U/N \rangle$', color='#ff7f0e', alpha=0.85, lw=1.3)
    ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=12, fontweight='bold')
    ax_bot.set_ylabel(r'Potential Energy $U / N$', fontsize=12, fontweight='bold')
    ax_bot.grid(True, linestyle=':', alpha=0.6)
    ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax_bot.set_ylim(-3.3, -2.4)
    
    # ------------------ COLUMN 3: SITE-SITE PATCHY (SSP) ------------------
    d_ssp = data['ssp']
    ax_top = axes[0, 2]
    ax_bot = axes[1, 2]
    
    # Top: Density
    ax_top.plot(d_ssp['d_npt'][:, 0], d_ssp['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
    ax_top.axhline(np.mean(d_ssp['d_nvt'][:, 2]), color='#d62728', linestyle='--', lw=2.0,
                   label=r'$NVT$ $\langle \rho^* \rangle = 0.7238$')
    ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=12, fontweight='bold')
    ax_top.set_title('(c) Site-Site Patchy (SSP)', fontsize=13, fontweight='bold', pad=8)
    ax_top.grid(True, linestyle=':', alpha=0.6)
    ax_top.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax_top.set_ylim(0.05, 0.78)
    
    # Bottom: Energy
    ax_bot.plot(d_ssp['d_npt'][:, 0], d_ssp['d_npt'][:, 7], label=r'$NpT$ $\langle U/N \rangle$', color='#2ca02c', alpha=0.85, lw=1.3)
    ax_bot.plot(d_ssp['d_nvt'][:, 0], d_ssp['d_nvt'][:, 7], label=r'$NVT$ $\langle U/N \rangle$', color='#ff7f0e', alpha=0.85, lw=1.3)
    ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=12, fontweight='bold')
    ax_bot.set_ylabel(r'Potential Energy $U / N$', fontsize=12, fontweight='bold')
    ax_bot.grid(True, linestyle=':', alpha=0.6)
    ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax_bot.set_ylim(-2.0, -1.2)
    
    out_pdf = os.path.join(out_dir, 'consistency_all.pdf')
    out_png = os.path.join(out_dir, 'consistency_all.png')
    plt.savefig(out_pdf, bbox_inches='tight')
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Grouped figure saved: {out_pdf} and {out_png}")

def generate_individual_figures(data, out_dir):
    # For each system, generate a 2x1 stacked figure sharing x-axis
    for sys_name in ['hs', 'lj', 'ssp']:
        fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(4.6, 5.8), sharex=True)
        plt.subplots_adjust(hspace=0.08, top=0.93, bottom=0.11, left=0.18, right=0.96)
        
        d = data[sys_name]
        
        if sys_name == 'hs':
            title = 'Hard Spheres (HS)'
            # Top: Density
            ax_top.plot(d['d_npt'][:, 0], d['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
            ax_top.axhline(np.mean(d['d_nvt'][:, 2]), color='#d62728', linestyle='--', lw=2.0,
                           label=r'$NVT$ $\langle \rho^* \rangle = 0.3364$')
            ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=11, fontweight='bold')
            ax_top.set_title(title, fontsize=12, fontweight='bold', pad=6)
            ax_top.grid(True, linestyle=':', alpha=0.6)
            ax_top.legend(loc='lower right', frameon=True, framealpha=0.9, fontsize=9)
            ax_top.set_ylim(0.295, 0.345)
            
            # Bottom: Virial Pressure
            s_npt, p_npt = d['s_npt'], d['p_npt']
            s_nvt, p_nvt = d['s_nvt'], d['p_nvt']
            eta_nvt = (np.pi / 6.0) * np.mean(d['d_nvt'][:, 2])
            p_cs = np.mean(d['d_nvt'][:, 2]) * (1.0 + eta_nvt + eta_nvt**2 - eta_nvt**3) / ((1.0 - eta_nvt)**3)
            
            if len(s_npt) > 0:
                ax_bot.plot(s_npt, p_npt, label=r'$NpT$ $P_{\rm vir}^*$', color='#2ca02c', alpha=0.75, lw=1.3)
            if len(s_nvt) > 0:
                ax_bot.plot(s_nvt, p_nvt, label=r'$NVT$ $P_{\rm vir}^*$', color='#ff7f0e', alpha=0.85, lw=1.3)
            ax_bot.axhline(p_cs, color='#d62728', linestyle='--', lw=2.0, label=f'CS EOS ({p_cs:.4f})')
            ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=11, fontweight='bold')
            ax_bot.set_ylabel(r'Virial Pressure $P\sigma^3/k_B T$', fontsize=11, fontweight='bold')
            ax_bot.grid(True, linestyle=':', alpha=0.6)
            ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9, fontsize=9)
            ax_bot.set_ylim(0.55, 0.80)
            
        elif sys_name == 'lj':
            title = 'Lennard-Jones (LJ)'
            # Top: Density
            ax_top.plot(d['d_npt'][:, 0], d['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
            ax_top.axhline(np.mean(d['d_nvt'][:, 2]), color='#d62728', linestyle='--', lw=2.0,
                           label=r'$NVT$ $\langle \rho^* \rangle = 0.4703$')
            ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=11, fontweight='bold')
            ax_top.set_title(title, fontsize=12, fontweight='bold', pad=6)
            ax_top.grid(True, linestyle=':', alpha=0.6)
            ax_top.legend(loc='upper right', frameon=True, framealpha=0.9, fontsize=9)
            ax_top.set_ylim(0.445, 0.585)
            
            # Bottom: Energy
            ax_bot.plot(d['d_npt'][:, 0], d['d_npt'][:, 7], label=r'$NpT$ $\langle U/N \rangle$', color='#2ca02c', alpha=0.85, lw=1.3)
            ax_bot.plot(d['d_nvt'][:, 0], d['d_nvt'][:, 7], label=r'$NVT$ $\langle U/N \rangle$', color='#ff7f0e', alpha=0.85, lw=1.3)
            ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=11, fontweight='bold')
            ax_bot.set_ylabel(r'Potential Energy $U / N$', fontsize=11, fontweight='bold')
            ax_bot.grid(True, linestyle=':', alpha=0.6)
            ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9, fontsize=9)
            ax_bot.set_ylim(-3.3, -2.4)
            
        elif sys_name == 'ssp':
            title = 'Site-Site Patchy (SSP)'
            # Top: Density
            ax_top.plot(d['d_npt'][:, 0], d['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
            ax_top.axhline(np.mean(d['d_nvt'][:, 2]), color='#d62728', linestyle='--', lw=2.0,
                           label=r'$NVT$ $\langle \rho^* \rangle = 0.7238$')
            ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=11, fontweight='bold')
            ax_top.set_title(title, fontsize=12, fontweight='bold', pad=6)
            ax_top.grid(True, linestyle=':', alpha=0.6)
            ax_top.legend(loc='lower right', frameon=True, framealpha=0.9, fontsize=9)
            ax_top.set_ylim(0.05, 0.78)
            
            # Bottom: Energy
            ax_bot.plot(d['d_npt'][:, 0], d['d_npt'][:, 7], label=r'$NpT$ $\langle U/N \rangle$', color='#2ca02c', alpha=0.85, lw=1.3)
            ax_bot.plot(d['d_nvt'][:, 0], d['d_nvt'][:, 7], label=r'$NVT$ $\langle U/N \rangle$', color='#ff7f0e', alpha=0.85, lw=1.3)
            ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=11, fontweight='bold')
            ax_bot.set_ylabel(r'Potential Energy $U / N$', fontsize=11, fontweight='bold')
            ax_bot.grid(True, linestyle=':', alpha=0.6)
            ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9, fontsize=9)
            ax_bot.set_ylim(-2.0, -1.2)
            
        out_pdf = os.path.join(out_dir, f'consistency_{sys_name}.pdf')
        out_png = os.path.join(out_dir, f'consistency_{sys_name}.png')
        plt.savefig(out_pdf, bbox_inches='tight')
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {out_pdf} and {out_png}")

if __name__ == '__main__':
    base_dir = '/home/e.lomba/MC_Checkerboard/test_cases/npt_nvt_consistency'
    out_dir = '/home/e.lomba/MC_Checkerboard/manuscript_cpc/figures'
    data = load_data(base_dir)
    generate_grouped_figure(data, out_dir)
    generate_individual_figures(data, out_dir)
