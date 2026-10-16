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
import sys
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

def load_data(base_dir):
    data = {}
    for sys_name in ['hs', 'lj', 'ssp']:
        sys_dir = os.path.join(base_dir, sys_name)
        npt_file = os.path.join(sys_dir, 'npt', 'run-data.dat')
        nvt_file = os.path.join(sys_dir, 'nvt', 'run-data.dat')
        d_npt = np.loadtxt(npt_file) if os.path.exists(npt_file) else np.empty((0, 8))
        d_nvt = np.loadtxt(nvt_file) if os.path.exists(nvt_file) else np.empty((0, 8))
        
        is_hs = (sys_name == 'hs')
        s_npt, p_npt = parse_virial_p(os.path.join(sys_dir, 'npt', 'output.out'), is_hs=is_hs)
        s_nvt, p_nvt = parse_virial_p(os.path.join(sys_dir, 'nvt', 'output.out'), is_hs=is_hs)

        data[sys_name] = {
            'd_npt': d_npt,
            'd_nvt': d_nvt,
            's_npt': s_npt,
            'p_npt': p_npt,
            's_nvt': s_nvt,
            'p_nvt': p_nvt,
        }
    return data

def generate_grouped_figure(data, out_dir):
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.2), sharex='col')
    plt.subplots_adjust(hspace=0.08, wspace=0.28, top=0.92, bottom=0.10, left=0.07, right=0.98)
    
    # ------------------ COLUMN 1: HARD SPHERES (HS) ------------------
    d_hs = data['hs']
    ax_top = axes[0, 0]
    ax_bot = axes[1, 0]
    
    # Top: Density
    rho_nvt_hs = np.mean(d_hs['d_nvt'][:, 2]) if len(d_hs['d_nvt']) > 0 else 0.3975
    ax_top.plot(d_hs['d_npt'][:, 0], d_hs['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
    ax_top.axhline(rho_nvt_hs, color='#d62728', linestyle='--', lw=2.0,
                   label=rf'$NVT$ $\langle \rho^* \rangle = {rho_nvt_hs:.4f}$')
    ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=12, fontweight='bold')
    ax_top.set_title('(a) Hard Spheres (HS)', fontsize=13, fontweight='bold', pad=8)
    ax_top.grid(True, linestyle=':', alpha=0.6)
    ax_top.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax_top.set_ylim(0.380, 0.415)
    
    # Bottom: Virial Pressure
    s_npt, p_npt = d_hs['s_npt'], d_hs['p_npt']
    s_nvt, p_nvt = d_hs['s_nvt'], d_hs['p_nvt']
    eta_nvt = (np.pi / 6.0) * rho_nvt_hs
    p_cs = rho_nvt_hs * (1.0 + eta_nvt + eta_nvt**2 - eta_nvt**3) / ((1.0 - eta_nvt)**3)
    
    if len(s_npt) > 0:
        ax_bot.plot(s_npt, p_npt, label=r'$NpT$ $P_{\rm vir}^*$', color='#2ca02c', alpha=0.75, lw=1.3)
    if len(s_nvt) > 0:
        ax_bot.plot(s_nvt, p_nvt, label=r'$NVT$ $P_{\rm vir}^*$', color='#ff7f0e', alpha=0.85, lw=1.3)
    ax_bot.axhline(p_cs, color='#d62728', linestyle='--', lw=2.0, label=f'CS EOS ({p_cs:.4f})')
    ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=12, fontweight='bold')
    ax_bot.set_ylabel(r'Virial Pressure $P\sigma^3/k_B T$', fontsize=12, fontweight='bold')
    ax_bot.grid(True, linestyle=':', alpha=0.6)
    ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax_bot.set_ylim(0.85, 1.15)
    
    # ------------------ COLUMN 2: LENNARD-JONES (LJ) ------------------
    d_lj = data['lj']
    ax_top = axes[0, 1]
    ax_bot = axes[1, 1]
    
    # Top: Density
    rho_nvt_lj = np.mean(d_lj['d_nvt'][:, 2]) if len(d_lj['d_nvt']) > 0 else 0.4785
    ax_top.plot(d_lj['d_npt'][:, 0], d_lj['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
    ax_top.axhline(rho_nvt_lj, color='#d62728', linestyle='--', lw=2.0,
                   label=rf'$NVT$ $\langle \rho^* \rangle = {rho_nvt_lj:.4f}$')
    ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=12, fontweight='bold')
    ax_top.set_title('(b) Lennard-Jones (LJ)', fontsize=13, fontweight='bold', pad=8)
    ax_top.grid(True, linestyle=':', alpha=0.6)
    ax_top.legend(loc='upper right', frameon=True, framealpha=0.9)
    ax_top.set_ylim(0.445, 0.610)
    
    # Bottom: Virial Pressure
    s_npt, p_npt = d_lj['s_npt'], d_lj['p_npt']
    s_nvt, p_nvt = d_lj['s_nvt'], d_lj['p_nvt']
    p_ref_lj = 0.60
    if len(s_npt) > 0:
        ax_bot.plot(s_npt, p_npt, label=r'$NpT$ $P_{\rm vir}^*$', color='#2ca02c', alpha=0.75, lw=1.3)
    if len(s_nvt) > 0:
        ax_bot.plot(s_nvt, p_nvt, label=r'$NVT$ $P_{\rm vir}^*$', color='#ff7f0e', alpha=0.85, lw=1.3)
    ax_bot.axhline(p_ref_lj, color='#d62728', linestyle='--', lw=2.0, label=f'Target $P^* = {p_ref_lj:.2f}$')
    ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=12, fontweight='bold')
    ax_bot.set_ylabel(r'Virial Pressure $P\sigma^3/\varepsilon$', fontsize=12, fontweight='bold')
    ax_bot.grid(True, linestyle=':', alpha=0.6)
    ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax_bot.set_ylim(0.45, 0.75)
    
    # ------------------ COLUMN 3: SITE-SITE PATCHY (SSP) ------------------
    d_ssp = data['ssp']
    ax_top = axes[0, 2]
    ax_bot = axes[1, 2]
    
    # Top: Density
    rho_nvt_ssp = np.mean(d_ssp['d_nvt'][:, 2]) if len(d_ssp['d_nvt']) > 0 else 0.7238
    ax_top.plot(d_ssp['d_npt'][:, 0], d_ssp['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
    ax_top.axhline(rho_nvt_ssp, color='#d62728', linestyle='--', lw=2.0,
                   label=rf'$NVT$ $\langle \rho^* \rangle = {rho_nvt_ssp:.4f}$')
    ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=12, fontweight='bold')
    ax_top.set_title('(c) Site-Site Patchy (SSP)', fontsize=13, fontweight='bold', pad=8)
    ax_top.grid(True, linestyle=':', alpha=0.6)
    ax_top.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax_top.set_ylim(0.745, 0.785)
    
    # Bottom: Virial Pressure
    s_npt, p_npt = d_ssp['s_npt'], d_ssp['p_npt']
    s_nvt, p_nvt = d_ssp['s_nvt'], d_ssp['p_nvt']
    p_ref_ssp = 0.05
    if len(s_npt) > 0:
        ax_bot.plot(s_npt, p_npt, label=r'$NpT$ $P_{\rm vir}^*$', color='#2ca02c', alpha=0.75, lw=1.3)
    if len(s_nvt) > 0:
        ax_bot.plot(s_nvt, p_nvt, label=r'$NVT$ $P_{\rm vir}^*$', color='#ff7f0e', alpha=0.85, lw=1.3)
    ax_bot.axhline(p_ref_ssp, color='#d62728', linestyle='--', lw=2.0, label=f'Target $P^* = {p_ref_ssp:.2f}$')
    ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=12, fontweight='bold')
    ax_bot.set_ylabel(r'Virial Pressure $P\sigma^3/\varepsilon$', fontsize=12, fontweight='bold')
    ax_bot.grid(True, linestyle=':', alpha=0.6)
    ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax_bot.set_ylim(-0.02, 0.12)
    
    out_pdf = os.path.join(out_dir, 'consistency_all.pdf')
    out_png = os.path.join(out_dir, 'consistency_all.png')
    plt.savefig(out_pdf, bbox_inches='tight')
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Grouped figure saved: {out_pdf} and {out_png}")

def generate_individual_figures(data, out_dir):
    for sys_name in ['hs', 'lj', 'ssp']:
        fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(4.6, 5.8), sharex=True)
        plt.subplots_adjust(hspace=0.08, top=0.93, bottom=0.11, left=0.18, right=0.96)
        
        d = data[sys_name]
        
        if sys_name == 'hs':
            title = 'Hard Spheres (HS)'
            rho_nvt_val = np.mean(d['d_nvt'][:, 2]) if len(d['d_nvt']) > 0 else 0.3975
            ax_top.plot(d['d_npt'][:, 0], d['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
            ax_top.axhline(rho_nvt_val, color='#d62728', linestyle='--', lw=2.0,
                           label=rf'$NVT$ $\langle \rho^* \rangle = {rho_nvt_val:.4f}$')
            ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=11, fontweight='bold')
            ax_top.set_title(title, fontsize=12, fontweight='bold', pad=6)
            ax_top.grid(True, linestyle=':', alpha=0.6)
            ax_top.legend(loc='lower right', frameon=True, framealpha=0.9, fontsize=9)
            ax_top.set_ylim(0.380, 0.415)
            
            s_npt, p_npt = d['s_npt'], d['p_npt']
            s_nvt, p_nvt = d['s_nvt'], d['p_nvt']
            eta_nvt = (np.pi / 6.0) * rho_nvt_val
            p_cs = rho_nvt_val * (1.0 + eta_nvt + eta_nvt**2 - eta_nvt**3) / ((1.0 - eta_nvt)**3)
            
            if len(s_npt) > 0:
                ax_bot.plot(s_npt, p_npt, label=r'$NpT$ $P_{\rm vir}^*$', color='#2ca02c', alpha=0.75, lw=1.3)
            if len(s_nvt) > 0:
                ax_bot.plot(s_nvt, p_nvt, label=r'$NVT$ $P_{\rm vir}^*$', color='#ff7f0e', alpha=0.85, lw=1.3)
            ax_bot.axhline(p_cs, color='#d62728', linestyle='--', lw=2.0, label=f'CS EOS ({p_cs:.4f})')
            ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=11, fontweight='bold')
            ax_bot.set_ylabel(r'Virial Pressure $P\sigma^3/k_B T$', fontsize=11, fontweight='bold')
            ax_bot.grid(True, linestyle=':', alpha=0.6)
            ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9, fontsize=9)
            ax_bot.set_ylim(0.85, 1.15)
            
        elif sys_name == 'lj':
            title = 'Lennard-Jones (LJ)'
            rho_nvt_val = np.mean(d['d_nvt'][:, 2]) if len(d['d_nvt']) > 0 else 0.4785
            ax_top.plot(d['d_npt'][:, 0], d['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
            ax_top.axhline(rho_nvt_val, color='#d62728', linestyle='--', lw=2.0,
                           label=rf'$NVT$ $\langle \rho^* \rangle = {rho_nvt_val:.4f}$')
            ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=11, fontweight='bold')
            ax_top.set_title(title, fontsize=12, fontweight='bold', pad=6)
            ax_top.grid(True, linestyle=':', alpha=0.6)
            ax_top.legend(loc='upper right', frameon=True, framealpha=0.9, fontsize=9)
            ax_top.set_ylim(0.445, 0.610)
            
            s_npt, p_npt = d['s_npt'], d['p_npt']
            s_nvt, p_nvt = d['s_nvt'], d['p_nvt']
            p_ref_lj = 0.60
            if len(s_npt) > 0:
                ax_bot.plot(s_npt, p_npt, label=r'$NpT$ $P_{\rm vir}^*$', color='#2ca02c', alpha=0.75, lw=1.3)
            if len(s_nvt) > 0:
                ax_bot.plot(s_nvt, p_nvt, label=r'$NVT$ $P_{\rm vir}^*$', color='#ff7f0e', alpha=0.85, lw=1.3)
            ax_bot.axhline(p_ref_lj, color='#d62728', linestyle='--', lw=2.0, label=f'Target $P^* = {p_ref_lj:.2f}$')
            ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=11, fontweight='bold')
            ax_bot.set_ylabel(r'Virial Pressure $P\sigma^3/\varepsilon$', fontsize=11, fontweight='bold')
            ax_bot.grid(True, linestyle=':', alpha=0.6)
            ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9, fontsize=9)
            ax_bot.set_ylim(0.45, 0.75)
            
        elif sys_name == 'ssp':
            title = 'Site-Site Patchy (SSP)'
            rho_nvt_val = np.mean(d['d_nvt'][:, 2]) if len(d['d_nvt']) > 0 else 0.7238
            ax_top.plot(d['d_npt'][:, 0], d['d_npt'][:, 2], label=r'$NpT$ $\rho^*(t)$', color='#1f77b4', lw=2.0)
            ax_top.axhline(rho_nvt_val, color='#d62728', linestyle='--', lw=2.0,
                           label=rf'$NVT$ $\langle \rho^* \rangle = {rho_nvt_val:.4f}$')
            ax_top.set_ylabel(r'Reduced Density $\rho^*$', fontsize=11, fontweight='bold')
            ax_top.set_title(title, fontsize=12, fontweight='bold', pad=6)
            ax_top.grid(True, linestyle=':', alpha=0.6)
            ax_top.legend(loc='lower right', frameon=True, framealpha=0.9, fontsize=9)
            ax_top.set_ylim(0.745, 0.785)
            
            s_npt, p_npt = d['s_npt'], d['p_npt']
            s_nvt, p_nvt = d['s_nvt'], d['p_nvt']
            p_ref_ssp = 0.05
            if len(s_npt) > 0:
                ax_bot.plot(s_npt, p_npt, label=r'$NpT$ $P_{\rm vir}^*$', color='#2ca02c', alpha=0.75, lw=1.3)
            if len(s_nvt) > 0:
                ax_bot.plot(s_nvt, p_nvt, label=r'$NVT$ $P_{\rm vir}^*$', color='#ff7f0e', alpha=0.85, lw=1.3)
            ax_bot.axhline(p_ref_ssp, color='#d62728', linestyle='--', lw=2.0, label=f'Target $P^* = {p_ref_ssp:.2f}$')
            ax_bot.set_xlabel('Monte Carlo Sweeps', fontsize=11, fontweight='bold')
            ax_bot.set_ylabel(r'Virial Pressure $P\sigma^3/\varepsilon$', fontsize=11, fontweight='bold')
            ax_bot.grid(True, linestyle=':', alpha=0.6)
            ax_bot.legend(loc='lower right', frameon=True, framealpha=0.9, fontsize=9)
            ax_bot.set_ylim(-0.02, 0.12)
            
        out_pdf = os.path.join(out_dir, f'consistency_{sys_name}.pdf')
        out_png = os.path.join(out_dir, f'consistency_{sys_name}.png')
        plt.savefig(out_pdf, bbox_inches='tight')
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {out_pdf} and {out_png}")

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.abspath(os.path.join(script_dir, '..', 'test_cases', 'npt_nvt_consistency'))
    out_dir = os.path.abspath(os.path.join(script_dir, 'figures'))
    if len(sys.argv) > 1:
        base_dir = sys.argv[1]
    if len(sys.argv) > 2:
        out_dir = sys.argv[2]
    data = load_data(base_dir)
    generate_grouped_figure(data, out_dir)
    generate_individual_figures(data, out_dir)
