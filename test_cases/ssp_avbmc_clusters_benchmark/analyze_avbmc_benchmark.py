#!/usr/bin/env python3
"""
==============================================================================
AVBMC vs Standard Canonical MC Comparative Benchmark Analysis Tool
==============================================================================
Parses simulation logs, cluster evolutions, and final cluster configurations
from:
  - run_with_avbmc/
  - run_without_avbmc/

Computes:
  1. Cluster size distribution (CSD) and mass distribution.
  2. Cluster growth/dissolution kinetics (N_clusters, S_max, % clustered).
  3. Acceptance rates (canonical translations/rotations vs AVBMC in/out).
  4. Computational throughput (ms/sweep, GPU execution time, speedup).
  5. Generates high-resolution multi-panel publication figures (PNG).

Author: Enrique Lomba / Antonio Diaz Pozuelo / Antigravity
Date: 2026
==============================================================================
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10.5,
    'ytick.labelsize': 10.5,
    'legend.fontsize': 10,
    'lines.linewidth': 2.0,
    'figure.dpi': 300
})

def parse_clusevol(filepath):
    """Parses clusevol_mc.dat file."""
    if not os.path.exists(filepath):
        return None
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    step = int(parts[0])
                    n_cl = int(parts[1])
                    max_s = int(parts[2])
                    n_clustered = int(parts[3])
                    pct = float(parts[4])
                    data.append([step, n_cl, max_s, n_clustered, pct])
                except ValueError:
                    continue
    return np.array(data) if data else None

def parse_lammpstrj_clusters(filepath):
    """
    Parses mclast_clconf.lammpstrj to extract exact cluster size distribution.
    Format: ITEM: ATOMS id mol type x y z ...
    where mol is the cluster ID.
    """
    if not os.path.exists(filepath):
        return None
    cluster_counts = {}
    reading_atoms = False
    with open(filepath, 'r') as f:
        for line in f:
            if 'ITEM: ATOMS' in line:
                reading_atoms = True
                continue
            if reading_atoms:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        mol_id = int(parts[1])
                        if mol_id > 0:
                            cluster_counts[mol_id] = cluster_counts.get(mol_id, 0) + 1
                    except ValueError:
                        continue
    sizes = list(cluster_counts.values())
    return np.array(sizes) if sizes else np.array([])

def parse_output_log(filepath):
    """Parses output.out to extract MC steps, energies, acceptance rates, and timings."""
    if not os.path.exists(filepath):
        return None
    
    steps = []
    energies = []
    p_trans = []
    p_rot = []
    p_av_in = []
    p_av_out = []
    cpu_times = []
    gpu_times = []

    total_cpu_time = None
    total_gpu_time = None

    with open(filepath, 'r') as f:
        for line in f:
            if 'Total CPU Time spent in main loop =' in line:
                try:
                    total_cpu_time = float(line.split('=')[1].replace('seconds.', '').strip())
                except:
                    pass
            elif 'Total GPU Time spent in main loop =' in line:
                try:
                    total_gpu_time = float(line.split('=')[1].replace('seconds.', '').strip())
                except:
                    pass
            
            # Line format with AVBMC:
            # Step En_tot/N En_tot/Ns P_Trans P_Rot N_Trans/step move/part N_Cell N_Clust Max_Cl %Clust P_AV_in P_AV_out CPU(s) GPU(s)
            parts = line.strip().split()
            if len(parts) >= 15:
                try:
                    step = int(parts[0])
                    en_n = float(parts[1])
                    pt = float(parts[3])
                    pr = float(parts[4])
                    pin = float(parts[11])
                    pout = float(parts[12])
                    cpu_t = float(parts[13])
                    gpu_t = float(parts[14])

                    steps.append(step)
                    energies.append(en_n)
                    p_trans.append(pt)
                    p_rot.append(pr)
                    p_av_in.append(pin)
                    p_av_out.append(pout)
                    cpu_times.append(cpu_t)
                    gpu_times.append(gpu_t)
                except ValueError:
                    continue
            elif len(parts) >= 13:
                # Standard format without AVBMC:
                try:
                    step = int(parts[0])
                    en_n = float(parts[1])
                    pt = float(parts[3])
                    pr = float(parts[4])
                    cpu_t = float(parts[11])
                    gpu_t = float(parts[12])

                    steps.append(step)
                    energies.append(en_n)
                    p_trans.append(pt)
                    p_rot.append(pr)
                    p_av_in.append(0.0)
                    p_av_out.append(0.0)
                    cpu_times.append(cpu_t)
                    gpu_times.append(gpu_t)
                except ValueError:
                    continue

    return {
        'step': np.array(steps),
        'energy_per_n': np.array(energies),
        'p_trans': np.array(p_trans),
        'p_rot': np.array(p_rot),
        'p_av_in': np.array(p_av_in),
        'p_av_out': np.array(p_av_out),
        'cpu_time': np.array(cpu_times),
        'gpu_time': np.array(gpu_times),
        'total_cpu_time': total_cpu_time,
        'total_gpu_time': total_gpu_time
    }

def run_analysis(base_dir='.'):
    dir_with = os.path.join(base_dir, 'run_with_avbmc')
    dir_without = os.path.join(base_dir, 'run_without_avbmc')

    print("=" * 80)
    print("  AVBMC VS CANONICAL MONTE CARLO BENCHMARK ANALYSIS")
    print("=" * 80)

    evol_with = parse_clusevol(os.path.join(dir_with, 'clusevol_mc.dat'))
    evol_without = parse_clusevol(os.path.join(dir_without, 'clusevol_mc.dat'))

    sizes_with = parse_lammpstrj_clusters(os.path.join(dir_with, 'mclast_clconf.lammpstrj'))
    sizes_without = parse_lammpstrj_clusters(os.path.join(dir_without, 'mclast_clconf.lammpstrj'))

    log_with = parse_output_log(os.path.join(dir_with, 'output.out'))
    log_without = parse_output_log(os.path.join(dir_without, 'output.out'))

    # Summary table output
    print(f"{'Metric':<35} | {'WITH AVBMC':<20} | {'WITHOUT AVBMC':<20}")
    print("-" * 80)

    if evol_with is not None and evol_without is not None:
        last_with = evol_with[-1]
        last_without = evol_without[-1]
        print(f"{'Final Cluster Count':<35} | {int(last_with[1]):<20} | {int(last_without[1]):<20}")
        print(f"{'Maximum Cluster Size (S_max)':<35} | {int(last_with[2]):<20} | {int(last_without[2]):<20}")
        print(f"{'Total Clustered Particles':<35} | {int(last_with[3]):<20} | {int(last_without[3]):<20}")
        print(f"{'Percentage Clustered (%)':<35} | {last_with[4]:<20.2f} | {last_without[4]:<20.2f}")

    if log_with and log_without:
        if len(log_with['p_av_in']) > 0:
            avg_p_in = np.mean(log_with['p_av_in'])
            avg_p_out = np.mean(log_with['p_av_out'])
            print(f"{'AVBMC In-Move Acceptance (P_in)':<35} | {avg_p_in*100:<19.2f}% | {'N/A':<20}")
            print(f"{'AVBMC Out-Move Acceptance (P_out)':<35} | {avg_p_out*100:<19.2f}% | {'N/A':<20}")

        if log_with['total_cpu_time'] and log_without['total_cpu_time']:
            print(f"{'Total CPU Time (s)':<35} | {log_with['total_cpu_time']:<20.2f} | {log_without['total_cpu_time']:<20.2f}")
        if log_with['total_gpu_time'] and log_without['total_gpu_time']:
            print(f"{'Total GPU Time (s)':<35} | {log_with['total_gpu_time']:<20.2f} | {log_without['total_gpu_time']:<20.2f}")
            
            n_steps_with = len(log_with['step']) * 5 if len(log_with['step']) > 0 else 200
            n_steps_without = len(log_without['step']) * 5 if len(log_without['step']) > 0 else 200
            ms_sweep_with = (log_with['total_gpu_time'] / n_steps_with) * 1000.0
            ms_sweep_without = (log_without['total_gpu_time'] / n_steps_without) * 1000.0
            print(f"{'GPU Time per MC Sweep (ms)':<35} | {ms_sweep_with:<20.2f} | {ms_sweep_without:<20.2f}")

    print("=" * 80)

    # ==============================================================================
    # 1. PRIMARY PUBLICATION FIGURE 4: Kinetics & Acceptance Rates (Stacked, Sharing X-axis)
    #    Eliminates old 4a (Cluster Abundance), stacks S_max (top) and Acceptance (bottom)
    # ==============================================================================
    fig_kin, (ax_smax, ax_acc) = plt.subplots(2, 1, figsize=(6.2, 7.6), sharex=True)
    plt.subplots_adjust(hspace=0.14)

    # Panel A: Maximum Cluster Size Evolution (formerly 4b)
    if evol_with is not None:
        ax_smax.plot(evol_with[:, 0], evol_with[:, 2], 'o-', color='#1f77b4', lw=2.2, label='With AVBMC')
    if evol_without is not None:
        ax_smax.plot(evol_without[:, 0], evol_without[:, 2], 's--', color='#d62728', lw=2.2, label='Without AVBMC (Canonical)')
    ax_smax.set_ylabel(r'Maximum Cluster Size $S_{\max}$', fontsize=12, fontweight='bold')
    ax_smax.set_title('(a) Maximum Cluster Growth Kinetics', fontsize=13, fontweight='bold')
    ax_smax.grid(True, linestyle=':', alpha=0.6)
    ax_smax.legend(fontsize=10.5, loc='center right')

    # Panel B: Monte Carlo Acceptance Rates (formerly 4d)
    if log_with and len(log_with['p_av_in']) > 0:
        ax_acc.plot(log_with['step'], log_with['p_av_in'] * 100, 'o-', color='#2ca02c', lw=2.0, label='AVBMC In-move Acc (%)')
        ax_acc.plot(log_with['step'], log_with['p_av_out'] * 100, 'v-', color='#ff7f0e', lw=2.0, label='AVBMC Out-move Acc (%)')
        ax_acc.plot(log_with['step'], log_with['p_trans'] * 100, 'k--', lw=1.5, alpha=0.7, label='Canonical Trans Acc (%)')
    ax_acc.set_xlabel('Monte Carlo Steps', fontsize=12, fontweight='bold')
    ax_acc.set_ylabel('Acceptance Rate (%)', fontsize=12, fontweight='bold')
    ax_acc.set_title('(b) Monte Carlo Acceptance Rates', fontsize=13, fontweight='bold')
    ax_acc.grid(True, linestyle=':', alpha=0.6)
    ax_acc.legend(fontsize=10.5, loc='center right')

    fig_kin.tight_layout()
    plot_kin_png = os.path.join(base_dir, 'avbmc_benchmark_comparison.png')
    plot_kin_pdf = os.path.join(base_dir, 'avbmc_benchmark_comparison.pdf')
    fig_kin.savefig(plot_kin_png, dpi=300, bbox_inches='tight')
    fig_kin.savefig(plot_kin_pdf, bbox_inches='tight')

    cpc_dir = os.path.abspath(os.path.join(base_dir, '..', '..', 'manuscript_cpc', 'figures'))
    if os.path.isdir(cpc_dir):
        fig_kin.savefig(os.path.join(cpc_dir, 'avbmc_benchmark_comparison.png'), dpi=300, bbox_inches='tight')
        fig_kin.savefig(os.path.join(cpc_dir, 'avbmc_benchmark_comparison.pdf'), bbox_inches='tight')
    plt.close(fig_kin)

    # ==============================================================================
    # 2. SEPARATED FIGURE: Final Cluster Size Distribution (CSD) (kept separated)
    # ==============================================================================
    fig_csd, ax_csd = plt.subplots(1, 1, figsize=(6.2, 4.6))
    if sizes_with is not None and len(sizes_with) > 0:
        bins = np.logspace(np.log10(max(1, np.min(sizes_with))), np.log10(np.max(sizes_with) + 1), 25)
        ax_csd.hist(sizes_with, bins=bins, color='#1f77b4', alpha=0.65, label='With AVBMC', edgecolor='black', density=True)
    if sizes_without is not None and len(sizes_without) > 0:
        bins = np.logspace(np.log10(max(1, np.min(sizes_without))), np.log10(np.max(sizes_without) + 1), 25)
        ax_csd.hist(sizes_without, bins=bins, color='#d62728', alpha=0.5, label='Without AVBMC', edgecolor='black', density=True)
    ax_csd.set_xscale('log')
    ax_csd.set_yscale('log')
    ax_csd.set_xlabel(r'Cluster Size $s$ (particles)', fontsize=12, fontweight='bold')
    ax_csd.set_ylabel(r'Probability Density $P(s)$', fontsize=12, fontweight='bold')
    ax_csd.set_title('Cluster Size Distribution (CSD)', fontsize=13, fontweight='bold')
    ax_csd.grid(True, linestyle=':', alpha=0.6)
    ax_csd.legend(fontsize=11)

    fig_csd.tight_layout()
    plot_csd_png = os.path.join(base_dir, 'avbmc_csd_distribution.png')
    plot_csd_pdf = os.path.join(base_dir, 'avbmc_csd_distribution.pdf')
    fig_csd.savefig(plot_csd_png, dpi=300, bbox_inches='tight')
    fig_csd.savefig(plot_csd_pdf, bbox_inches='tight')
    if os.path.isdir(cpc_dir):
        fig_csd.savefig(os.path.join(cpc_dir, 'avbmc_csd_distribution.png'), dpi=300, bbox_inches='tight')
        fig_csd.savefig(os.path.join(cpc_dir, 'avbmc_csd_distribution.pdf'), bbox_inches='tight')
    plt.close(fig_csd)

    # ==============================================================================
    # 3. UNIFIED 3-PANEL FIGURE: Left = stacked a & b sharing x-axis, Right = c (CSD) separated
    # ==============================================================================
    fig_all = plt.figure(figsize=(12, 6.2))
    gs = fig_all.add_gridspec(2, 2, width_ratios=[1.0, 1.0], hspace=0.18, wspace=0.25)
    ax1_all = fig_all.add_subplot(gs[0, 0])
    ax2_all = fig_all.add_subplot(gs[1, 0], sharex=ax1_all)
    ax3_all = fig_all.add_subplot(gs[:, 1])

    # Subplot A: S_max
    if evol_with is not None:
        ax1_all.plot(evol_with[:, 0], evol_with[:, 2], 'o-', color='#1f77b4', lw=2.2, label='With AVBMC')
    if evol_without is not None:
        ax1_all.plot(evol_without[:, 0], evol_without[:, 2], 's--', color='#d62728', lw=2.2, label='Without AVBMC')
    ax1_all.set_ylabel(r'Max Cluster Size $S_{\max}$', fontsize=11, fontweight='bold')
    ax1_all.set_title('(a) Cluster Growth Kinetics', fontsize=12, fontweight='bold')
    ax1_all.grid(True, linestyle=':', alpha=0.6)
    ax1_all.legend(fontsize=9.5, loc='center right')
    plt.setp(ax1_all.get_xticklabels(), visible=False)

    # Subplot B: Acceptance
    if log_with and len(log_with['p_av_in']) > 0:
        ax2_all.plot(log_with['step'], log_with['p_av_in'] * 100, 'o-', color='#2ca02c', lw=2.0, label='AVBMC In (%)')
        ax2_all.plot(log_with['step'], log_with['p_av_out'] * 100, 'v-', color='#ff7f0e', lw=2.0, label='AVBMC Out (%)')
        ax2_all.plot(log_with['step'], log_with['p_trans'] * 100, 'k--', lw=1.5, alpha=0.7, label='Canonical Trans (%)')
    ax2_all.set_xlabel('Monte Carlo Steps', fontsize=11, fontweight='bold')
    ax2_all.set_ylabel('Acceptance (%)', fontsize=11, fontweight='bold')
    ax2_all.set_title('(b) Acceptance Rates', fontsize=12, fontweight='bold')
    ax2_all.grid(True, linestyle=':', alpha=0.6)
    ax2_all.legend(fontsize=9.5, loc='center right')

    # Subplot C: CSD
    if sizes_with is not None and len(sizes_with) > 0:
        bins = np.logspace(np.log10(max(1, np.min(sizes_with))), np.log10(np.max(sizes_with) + 1), 25)
        ax3_all.hist(sizes_with, bins=bins, color='#1f77b4', alpha=0.65, label='With AVBMC', edgecolor='black', density=True)
    if sizes_without is not None and len(sizes_without) > 0:
        bins = np.logspace(np.log10(max(1, np.min(sizes_without))), np.log10(np.max(sizes_without) + 1), 25)
        ax3_all.hist(sizes_without, bins=bins, color='#d62728', alpha=0.5, label='Without AVBMC', edgecolor='black', density=True)
    ax3_all.set_xscale('log')
    ax3_all.set_yscale('log')
    ax3_all.set_xlabel(r'Cluster Size $s$ (particles)', fontsize=11, fontweight='bold')
    ax3_all.set_ylabel(r'Probability Density $P(s)$', fontsize=11, fontweight='bold')
    ax3_all.set_title('(c) Cluster Size Distribution (CSD)', fontsize=12, fontweight='bold')
    ax3_all.grid(True, linestyle=':', alpha=0.6)
    ax3_all.legend(fontsize=10.5)

    plot_all_png = os.path.join(base_dir, 'avbmc_benchmark_all.png')
    plot_all_pdf = os.path.join(base_dir, 'avbmc_benchmark_all.pdf')
    fig_all.savefig(plot_all_png, dpi=300, bbox_inches='tight')
    fig_all.savefig(plot_all_pdf, bbox_inches='tight')
    if os.path.isdir(cpc_dir):
        fig_all.savefig(os.path.join(cpc_dir, 'avbmc_benchmark_all.png'), dpi=300, bbox_inches='tight')
        fig_all.savefig(os.path.join(cpc_dir, 'avbmc_benchmark_all.pdf'), bbox_inches='tight')
    plt.close(fig_all)
    print(f"  >>> Benchmark comparison plots saved: {plot_kin_png}, {plot_csd_png}, and {plot_all_png}")

if __name__ == '__main__':
    run_analysis(base_dir='.')
