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

    # Generate Publication Figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    plt.subplots_adjust(hspace=0.32, wspace=0.25)

    # 1. Cluster Count Evolution
    ax1 = axes[0, 0]
    if evol_with is not None:
        ax1.plot(evol_with[:, 0], evol_with[:, 1], 'o-', color='#1f77b4', lw=2.2, label='With AVBMC')
    if evol_without is not None:
        ax1.plot(evol_without[:, 0], evol_without[:, 1], 's--', color='#d62728', lw=2.2, label='Without AVBMC (Canonical)')
    ax1.set_xlabel('Monte Carlo Steps', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Clusters', fontsize=12, fontweight='bold')
    ax1.set_title('(a) Cluster Abundance Kinetics', fontsize=13, fontweight='bold')
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(fontsize=11)

    # 2. Maximum Cluster Size Evolution
    ax2 = axes[0, 1]
    if evol_with is not None:
        ax2.plot(evol_with[:, 0], evol_with[:, 2], 'o-', color='#1f77b4', lw=2.2, label='With AVBMC')
    if evol_without is not None:
        ax2.plot(evol_without[:, 0], evol_without[:, 2], 's--', color='#d62728', lw=2.2, label='Without AVBMC (Canonical)')
    ax2.set_xlabel('Monte Carlo Steps', fontsize=12, fontweight='bold')
    ax2.set_ylabel(r'Maximum Cluster Size $S_{\max}$', fontsize=12, fontweight='bold')
    ax2.set_title('(b) Maximum Cluster Growth Kinetics', fontsize=13, fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(fontsize=11)

    # 3. Cluster Size Distribution (Histogram)
    ax3 = axes[1, 0]
    if sizes_with is not None and len(sizes_with) > 0:
        bins = np.logspace(np.log10(max(1, np.min(sizes_with))), np.log10(np.max(sizes_with) + 1), 25)
        ax3.hist(sizes_with, bins=bins, color='#1f77b4', alpha=0.65, label='With AVBMC', edgecolor='black', density=True)
    if sizes_without is not None and len(sizes_without) > 0:
        bins = np.logspace(np.log10(max(1, np.min(sizes_without))), np.log10(np.max(sizes_without) + 1), 25)
        ax3.hist(sizes_without, bins=bins, color='#d62728', alpha=0.5, label='Without AVBMC', edgecolor='black', density=True)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel(r'Cluster Size $s$ (particles)', fontsize=12, fontweight='bold')
    ax3.set_ylabel(r'Probability Density $P(s)$', fontsize=12, fontweight='bold')
    ax3.set_title('(c) Cluster Size Distribution (CSD)', fontsize=13, fontweight='bold')
    ax3.grid(True, linestyle=':', alpha=0.6)
    ax3.legend(fontsize=11)

    # 4. Acceptance Rates & Performance
    ax4 = axes[1, 1]
    if log_with and len(log_with['p_av_in']) > 0:
        ax4.plot(log_with['step'], log_with['p_av_in'] * 100, 'o-', color='#2ca02c', lw=2.0, label='AVBMC In-move Acc (%)')
        ax4.plot(log_with['step'], log_with['p_av_out'] * 100, 'v-', color='#ff7f0e', lw=2.0, label='AVBMC Out-move Acc (%)')
        ax4.plot(log_with['step'], log_with['p_trans'] * 100, 'k--', lw=1.5, alpha=0.7, label='Canonical Trans Acc (%)')
    ax4.set_xlabel('Monte Carlo Steps', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Acceptance Rate (%)', fontsize=12, fontweight='bold')
    ax4.set_title('(d) Monte Carlo Acceptance Rates', fontsize=13, fontweight='bold')
    ax4.grid(True, linestyle=':', alpha=0.6)
    ax4.legend(fontsize=11)

    plot_file = os.path.join(base_dir, 'avbmc_benchmark_comparison.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  >>> Benchmark comparison plot saved: {plot_file}")

if __name__ == '__main__':
    run_analysis(base_dir='.')
