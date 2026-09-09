#!/usr/bin/env python3
"""
plot_scaling_benchmark.py
Generates publication-quality dual-panel scaling and speedup figure for the CPC manuscript.
"""

import matplotlib.pyplot as plt
import numpy as np

# Data from Table 2 of benchmark.tex
N = np.array([8000, 15625, 32768, 125000])
mu_GPU = np.array([3.375, 1.728, 1.953, 2.744])

# CPU times (seconds)
t_cpu_nvt = np.array([9.95, 23.59, 38.63, 175.91])
t_cpu_npt = np.array([15.16, 35.95, 59.23, 254.09])

# GPU raw times (seconds)
t_gpu_raw_nvt = np.array([4.97, 8.73, 8.05, 7.86])
t_gpu_raw_npt = np.array([6.26, 11.88, 9.94, 10.14])

# GPU normalized times (seconds)
t_gpu_norm_nvt = t_gpu_raw_nvt / mu_GPU
t_gpu_norm_npt = t_gpu_raw_npt / mu_GPU

# Raw speedups
s_raw_nvt = t_cpu_nvt / t_gpu_raw_nvt
s_raw_npt = t_cpu_npt / t_gpu_raw_npt

# Effective throughput speedups
s_eff_nvt = t_cpu_nvt / t_gpu_norm_nvt
s_eff_npt = t_cpu_npt / t_gpu_norm_npt

# Style configuration
plt.style.use('seaborn-v0_8-paper' if 'seaborn-v0_8-paper' in plt.style.available else 'default')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8), dpi=300)

colors = {
    'cpu_nvt': '#1f77b4',
    'cpu_npt': '#0d47a1',
    'gpu_raw_nvt': '#2ca02c',
    'gpu_raw_npt': '#1b5e20',
    'gpu_norm_nvt': '#ff7f0e',
    'gpu_norm_npt': '#d84315',
    'eff_nvt': '#d62728',
    'eff_npt': '#b71c1c',
}

# ----------------- PANEL A: Execution Times -----------------
ax1.loglog(N, t_cpu_nvt, 's--', color=colors['cpu_nvt'], label=r'gpMC CPU ($NVT$)', linewidth=1.8, markersize=7)
ax1.loglog(N, t_cpu_npt, 'o--', color=colors['cpu_npt'], label=r'gpMC CPU ($NpT$)', linewidth=1.8, markersize=7)
ax1.loglog(N, t_gpu_raw_nvt, '^-.', color=colors['gpu_raw_nvt'], label=r'MCCB-gpu Raw ($NVT$)', linewidth=1.8, markersize=7)
ax1.loglog(N, t_gpu_raw_npt, 'v-.', color=colors['gpu_raw_npt'], label=r'MCCB-gpu Raw ($NpT$)', linewidth=1.8, markersize=7)
ax1.loglog(N, t_gpu_norm_nvt, 'D-', color=colors['gpu_norm_nvt'], label=r'MCCB-gpu Norm ($NVT$)', linewidth=2.0, markersize=7)
ax1.loglog(N, t_gpu_norm_npt, 'h-', color=colors['gpu_norm_npt'], label=r'MCCB-gpu Norm ($NpT$)', linewidth=2.0, markersize=7)

ax1.set_xlabel(r'Particle Count ($N$)', fontsize=12, fontweight='bold')
ax1.set_ylabel(r'Execution Time for 100 Sweeps (s)', fontsize=12, fontweight='bold')
ax1.set_title('(a) Execution Time Scaling', fontsize=13, fontweight='bold')
ax1.set_xticks(N)
ax1.set_xticklabels([r'8k', r'15.6k', r'32.8k', r'125k'])
ax1.grid(True, which='both', linestyle=':', alpha=0.6)
ax1.legend(loc='upper left', frameon=True, fontsize=8.5, ncol=2)

# ----------------- PANEL B: Speedup -----------------
ax2.plot(N, s_eff_nvt, 'D-', color=colors['eff_nvt'], label=r'Effective Throughput ($NVT$)', linewidth=2.2, markersize=8)
ax2.plot(N, s_eff_npt, 'h-', color=colors['eff_npt'], label=r'Effective Throughput ($NpT$)', linewidth=2.2, markersize=8)
ax2.plot(N, s_raw_nvt, '^--', color=colors['gpu_raw_nvt'], label=r'Raw Cycle ($NVT$)', linewidth=1.8, markersize=7)
ax2.plot(N, s_raw_npt, 'v--', color=colors['gpu_raw_npt'], label=r'Raw Cycle ($NpT$)', linewidth=1.8, markersize=7)

ax2.set_xscale('log')
ax2.set_xlabel(r'Particle Count ($N$)', fontsize=12, fontweight='bold')
ax2.set_ylabel(r'Speedup Factor ($S = t_{\mathrm{CPU}} / t_{\mathrm{GPU}}$)', fontsize=12, fontweight='bold')
ax2.set_title('(b) Speedup Factor across Ensembles', fontsize=13, fontweight='bold')
ax2.set_xticks(N)
ax2.set_xticklabels([r'8k', r'15.6k', r'32.8k', r'125k'])
ax2.grid(True, which='both', linestyle=':', alpha=0.6)
ax2.legend(loc='upper left', frameon=True, fontsize=9)

# Annotate peak speedup
ax2.annotate(r'$\mathbf{68.8\times}$ ($NpT$)',
             xy=(125000, s_eff_npt[-1]), xytext=(70000, 72),
             arrowprops=dict(facecolor='black', shrink=0.08, width=1, headwidth=6),
             fontsize=10, fontweight='bold', color=colors['eff_npt'])

ax2.annotate(r'$\mathbf{61.4\times}$ ($NVT$)',
             xy=(125000, s_eff_nvt[-1]), xytext=(70000, 50),
             arrowprops=dict(facecolor='black', shrink=0.08, width=1, headwidth=6),
             fontsize=10, fontweight='bold', color=colors['eff_nvt'])

plt.tight_layout()
out_path = '/home/e.lomba/MC_Checkerboard/manuscript_cpc/figures/scaling_benchmark.png'
plt.savefig(out_path, dpi=300)
print(f"Figure saved successfully to {out_path}")
