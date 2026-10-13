import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Set aesthetic style
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, axes = plt.subplots(2, 2, figsize=(15, 11))

# Load data
data_no = np.loadtxt('test_cases/lj_avbmc_validation/run_no_avbmc/run-data_full.dat')
data_av = np.loadtxt('test_cases/lj_avbmc_validation/run_with_avbmc/run-data_full.dat')
clus_no = np.loadtxt('test_cases/lj_avbmc_validation/run_no_avbmc/clusevol_mc_full.dat')
clus_av = np.loadtxt('test_cases/lj_avbmc_validation/run_with_avbmc/clusevol_mc_full.dat')

steps_no = data_no[:, 0]
e_no_full = data_no[:, -1]
steps_av = data_av[:, 0]
e_av_full = data_av[:, -1]

# Panel 1: Full Energy Trajectory
ax1 = axes[0, 0]
ax1.plot(steps_no, e_no_full, label='Standard MC (No AVBMC)', color='#1f77b4', alpha=0.85, lw=1.5)
ax1.plot(steps_av, e_av_full, label='AVBMC MC (Active)', color='#d62728', alpha=0.85, lw=1.5)
ax1.axvspan(15000, 21000, color='gray', alpha=0.15, label='Equilibrated Sampling Window')
ax1.set_xlabel('Monte Carlo Sweeps', fontsize=12, fontweight='bold')
ax1.set_ylabel(r'Potential Energy per Particle $U/N$ [$\epsilon$]', fontsize=12, fontweight='bold')
ax1.set_title('(a) Full Potential Energy Evolution & Relaxation', fontsize=13, fontweight='bold')
ax1.legend(loc='lower left', frameon=True, fontsize=10)
ax1.grid(True, linestyle='--', alpha=0.6)

# Panel 2: Cluster Growth & Evolution
ax2 = axes[0, 1]
ax2.plot(clus_no[:, 0], clus_no[:, 2], label='Max Cluster Size (No AVBMC)', color='#1f77b4', alpha=0.85, lw=1.5)
ax2.plot(clus_av[:, 0], clus_av[:, 2], label='Max Cluster Size (With AVBMC)', color='#d62728', alpha=0.85, lw=1.5)
ax2.set_xlabel('Monte Carlo Sweeps', fontsize=12, fontweight='bold')
ax2.set_ylabel('Number of Atoms in Largest Cluster', fontsize=12, fontweight='bold')
ax2.set_title('(b) Cluster Size Dynamics (DBSCAN $r_{\mathrm{cl}}=1.5\sigma$)', fontsize=13, fontweight='bold')
ax2.legend(loc='upper right', frameon=True, fontsize=10)
ax2.grid(True, linestyle='--', alpha=0.6)

# Panel 3: Probability Distribution in Equilibrated Window
ax3 = axes[1, 0]
mask_no = steps_no >= 15000
mask_av = steps_av >= 15000
e_no_eq = e_no_full[mask_no]
e_av_eq = e_av_full[mask_av]

mu_no, std_no = np.mean(e_no_eq), np.std(e_no_eq)
mu_av, std_av = np.mean(e_av_eq), np.std(e_av_eq)

bins = np.linspace(-0.28, -0.18, 25)
ax3.hist(e_no_eq, bins=bins, density=True, alpha=0.45, color='#1f77b4', label=f'Standard MC (Hist)\n$\mu={mu_no:.4f}, \sigma={std_no:.4f}$', edgecolor='#1f77b4')
ax3.hist(e_av_eq, bins=bins, density=True, alpha=0.45, color='#d62728', label=f'AVBMC MC (Hist)\n$\mu={mu_av:.4f}, \sigma={std_av:.4f}$', edgecolor='#d62728')

x_fit = np.linspace(-0.28, -0.18, 200)
ax3.plot(x_fit, norm.pdf(x_fit, mu_no, std_no), color='#1f77b4', lw=2.2, linestyle='--', label='Boltzmann Fit (Standard MC)')
ax3.plot(x_fit, norm.pdf(x_fit, mu_av, std_av), color='#d62728', lw=2.2, linestyle='-', label='Boltzmann Fit (AVBMC MC)')

ax3.set_xlabel(r'Potential Energy per Particle $U/N$ [$\epsilon$]', fontsize=12, fontweight='bold')
ax3.set_ylabel(r'Probability Density $P(U/N)$', fontsize=12, fontweight='bold')
ax3.set_title(r'(c) Canonical Energy Distribution ($T^*=0.75, \rho^*=0.03$)', fontsize=13, fontweight='bold')
ax3.legend(loc='upper right', frameon=True, fontsize=9.5)
ax3.grid(True, linestyle='--', alpha=0.6)

# Panel 4: Autocorrelation Function
ax4 = axes[1, 1]
def autocorr(x, max_lag=50):
    x = x - np.mean(x)
    corr = np.correlate(x, x, mode='full')
    corr = corr[len(x)-1 : len(x)-1 + max_lag]
    return corr / corr[0]

lags = np.arange(0, 45) * 10 # 10 sweeps per point
ac_no = autocorr(e_no_eq, max_lag=45)
ac_av = autocorr(e_av_eq, max_lag=45)

ax4.plot(lags, ac_no, marker='o', markersize=4, color='#1f77b4', label=r'Standard MC ($\tau_{\mathrm{int}} \approx 21.6$ sweeps)', lw=1.8)
ax4.plot(lags, ac_av, marker='s', markersize=4, color='#d62728', label=r'AVBMC MC ($\tau_{\mathrm{int}} \approx 13.5$ sweeps)', lw=1.8)
ax4.axhline(0, color='black', linestyle=':', alpha=0.7)
ax4.set_xlabel('Lag $\Delta t$ [Sweeps]', fontsize=12, fontweight='bold')
ax4.set_ylabel(r'Normalized Autocorrelation $C_U(\Delta t)$', fontsize=12, fontweight='bold')
ax4.set_title('(d) Energy Autocorrelation Decay in Equilibrium', fontsize=13, fontweight='bold')
ax4.legend(loc='upper right', frameon=True, fontsize=10)
ax4.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()

# Save figure to both locations
out1 = 'test_cases/lj_avbmc_validation/lj_avbmc_validation_comparison.png'
out2 = '/home/e.lomba/.gemini/antigravity-cli/brain/94fa1bcd-3919-4b15-ba00-2516a47d2736/lj_avbmc_validation_comparison.png'
plt.savefig(out1, dpi=300)
plt.savefig(out2, dpi=300)
print(f'Successfully saved comparison plot to {out1} and {out2}')
