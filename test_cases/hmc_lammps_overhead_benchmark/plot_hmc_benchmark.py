#!/usr/bin/env python3
"""
Publication-Quality Analysis & Plotting Script for HMC Overhead Benchmark
Generates a 4-panel figure for JCTC manuscript submission (No pandas dependency).
"""

import matplotlib.pyplot as plt
import numpy as np
import csv
import os

# Set publication style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'lines.linewidth': 2.0,
    'lines.markersize': 7,
    'figure.dpi': 300
})

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# 1. Load CSV data using standard csv module
def load_csv_data(filepath):
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        header = [h.strip() for h in next(reader)]
        rows = [[val.strip() for val in r] for r in reader if r]
    return header, rows

header_trans, rows_trans = load_csv_data(os.path.join(SCRIPT_DIR, 'data_transfer_breakdown.csv'))
header_ana, rows_ana = load_csv_data(os.path.join(SCRIPT_DIR, 'trajectory_scaling_analytic.csv'))
header_tab, rows_tab = load_csv_data(os.path.join(SCRIPT_DIR, 'trajectory_scaling_table.csv'))

# Parse numerical columns
nmd_arr = np.array([float(r[0]) for r in rows_ana])
time_pure_ana = np.array([float(r[1]) for r in rows_ana])
time_pers_ana = np.array([float(r[2]) for r in rows_ana])
time_base_ana = np.array([float(r[3]) for r in rows_ana])
ovhd_base_ana = np.array([float(r[4]) for r in rows_ana])
ovhd_pers_ana = np.array([float(r[5]) for r in rows_ana])

time_pure_tab = np.array([float(r[1]) for r in rows_tab])
time_pers_tab = np.array([float(r[2]) for r in rows_tab])
time_base_tab = np.array([float(r[3]) for r in rows_tab])
ovhd_base_tab = np.array([float(r[4]) for r in rows_tab])
ovhd_pers_tab = np.array([float(r[5]) for r in rows_tab])

# Transfer times
trans_dict = {r[0]: float(r[2]) for r in rows_trans}

# Create 2x2 subplot figure
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# ------------------------------------------------------------------------------
# PANEL A: Data Transfer & Conversion Overhead Breakdown (ms)
# ------------------------------------------------------------------------------
ax = axes[0, 0]

categories = [
    'GPU $\\leftrightarrow$ Host Sync\n(Convert CB / Upload)',
    'Rigid Body $\\leftrightarrow$ 5N Sites\n(Rotations / Triad)',
    'Disk File I/O\n(Write + Read data)',
    'In-Memory API\n(Scatter + Gather)'
]

times = [
    trans_dict.get('T1_CB_to_Host', 0.102) + trans_dict.get('T8_Host_to_CB', 3.419),
    trans_dict.get('T2_Patch_Gen', 0.190) + trans_dict.get('T7_Patch_Recon', 0.509),
    trans_dict.get('T3_Disk_Write', 47.865) + trans_dict.get('T4_Disk_Read', 75.010),
    trans_dict.get('T5_Scatter', 0.164) + trans_dict.get('T6_Gather', 0.264)
]

colors = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd']
bars = ax.barh(categories, times, color=colors, edgecolor='black', height=0.55, alpha=0.85)

# Annotate bars
for bar, t in zip(bars, times):
    width = bar.get_width()
    ax.text(width + 1.5, bar.get_y() + bar.get_height()/2.0, f'{t:.2f} ms',
            ha='left', va='center', fontweight='bold', fontsize=10)

ax.set_xlim(0, max(times) * 1.25)
ax.set_xlabel('Time per HMC Trial (ms)')
ax.set_title('(a) Data Transfer & Conversion Layer Breakdown', pad=10)
ax.grid(True, linestyle='--', alpha=0.4, axis='x')
ax.annotate('287x Speedup\n(RAM vs Disk)', xy=(times[3], 3), xytext=(25, 2.7),
            arrowprops=dict(facecolor='black', shrink=0.08, width=1, headwidth=6),
            fontsize=10, fontweight='bold', color='#2e7d32')

# ------------------------------------------------------------------------------
# PANEL B: Execution Time vs Trajectory Length (Analytic vs Tabulated)
# ------------------------------------------------------------------------------
ax = axes[0, 1]

ax.plot(nmd_arr, time_pure_ana, 'o--', color='#d62728', label='Analytic: Pure MD')
ax.plot(nmd_arr, time_base_ana, 's-', color='#b71c1c', label='Analytic: Reset (Baseline)')
ax.plot(nmd_arr, time_pure_tab, '^--', color='#1f77b4', label='Tabulated: Pure MD')
ax.plot(nmd_arr, time_base_tab, 'd-', color='#0d47a1', label='Tabulated: Reset (Baseline)')

ax.set_xlabel('MD Steps per HMC Trial ($N_{\\mathrm{md}}$)')
ax.set_ylabel('Execution Time per Step (s)')
ax.set_title('(b) Trajectory Execution Scaling (RTX PRO 4500)', pad=10)
ax.legend(frameon=True, facecolor='white', framealpha=0.9, loc='upper left')
ax.grid(True, linestyle='--', alpha=0.4)

# Speedup annotation
speedup_50 = time_pure_ana[3] / time_pure_tab[3]
speedup_500 = time_pure_ana[6] / time_pure_tab[6]
ax.text(230, 15, f'GPU Table Speedup:\n{speedup_50:.1f}x at $N_{{md}}=50$\n{speedup_500:.1f}x at $N_{{md}}=500$',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#e3f2fd', edgecolor='#1e88e5', alpha=0.9),
        fontsize=10, fontweight='bold')

# ------------------------------------------------------------------------------
# PANEL C: Overhead Percentage vs Trajectory Length
# ------------------------------------------------------------------------------
ax = axes[1, 0]

ax.plot(nmd_arr[1:], ovhd_base_tab[1:], 'o-', color='#e65100', label='Tabulated (Reset / File I/O)')
ax.plot(nmd_arr[1:], ovhd_base_ana[1:], 's-', color='#c2185b', label='Analytic (Reset / File I/O)')
ax.plot(nmd_arr[1:], ovhd_pers_tab[1:], '^--', color='#2e7d32', label='Tabulated (Persistent / In-Memory)')
ax.plot(nmd_arr[1:], ovhd_pers_ana[1:], 'v--', color='#00838f', label='Analytic (Persistent / In-Memory)')

# 50% crossover line
ax.axhline(50, color='gray', linestyle=':', alpha=0.7)
ax.text(320, 52, '50% Overhead Threshold', color='gray', fontsize=9, style='italic')

ax.set_xlabel('MD Steps per HMC Trial ($N_{\\mathrm{md}}$)')
ax.set_ylabel('Overhead Fraction (%)')
ax.set_title('(c) Task & Transfer Overhead vs MD Trajectory Length', pad=10)
ax.set_ylim(0, 85)
ax.legend(frameon=True, facecolor='white', framealpha=0.9, loc='upper right')
ax.grid(True, linestyle='--', alpha=0.4)

# ------------------------------------------------------------------------------
# PANEL D: Task Invocation Paradigm Comparison (Nmd = 50)
# ------------------------------------------------------------------------------
ax = axes[1, 1]

paradigms = ['Process Fork\n(External lmp)', 'Library Reset\n(clear + read_data)', 'Library Persistent\n(scatter_atoms)']
t_analytic_50 = [1.583, time_base_ana[3], time_pers_ana[3]]
t_table_50 = [1.198, time_base_tab[3], time_pers_tab[3]]

x = np.arange(len(paradigms))
width = 0.35

rects1 = ax.bar(x - width/2, t_analytic_50, width, label='Analytic Potential', color='#e57373', edgecolor='black', alpha=0.9)
rects2 = ax.bar(x + width/2, t_table_50, width, label='Tabulated Potential', color='#64b5f6', edgecolor='black', alpha=0.9)

ax.set_ylabel('Wall-Clock Time per Step (s)')
ax.set_title('(d) Task Invocation Mechanisms ($N_{\\mathrm{md}} = 50$)', pad=10)
ax.set_xticks(x)
ax.set_xticklabels(paradigms)
ax.legend(frameon=True, facecolor='white', framealpha=0.9)
ax.grid(True, linestyle='--', alpha=0.4, axis='y')

# Value labels on bars
for rect in rects1:
    h = rect.get_height()
    ax.annotate(f'{h:.2f}s', xy=(rect.get_x() + rect.get_width() / 2, h),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9, fontweight='bold')
for rect in rects2:
    h = rect.get_height()
    ax.annotate(f'{h:.2f}s', xy=(rect.get_x() + rect.get_width() / 2, h),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.tight_layout()
output_fig = os.path.join(SCRIPT_DIR, 'hmc_overhead_benchmark_comparison.png')
plt.savefig(output_fig, dpi=300)
print(f"Figure saved successfully to: {output_fig}")
