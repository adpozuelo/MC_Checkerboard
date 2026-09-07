# AVBMC Patchy Colloidal Clusters Benchmark: Detailed Balance & Performance

## Overview

This benchmark directory contains a standalone test case and verification suite for the **Association Volume Bias Monte Carlo (AVBMC)** implementation in `MC_Checkerboard` for cluster-forming systems, prepared for submission to the *Journal of Chemical Theory and Computation* (JCTC).

The physical test case is based on `examples/SSP`: a site-site patchy colloidal fluid (Palaia model) with four tetrahedral patches, diluted by a factor of 100 ($\rho \to \rho/100 = 0.001218$, $N = 7695$ particles, $L = 184.8576\sigma$) at $T^* = 0.1540$. Under these conditions, the system aggregates into finite-size colloidal clusters surrounded by an ultra-dilute gas phase.

---

## Directory Structure

- `data_box_expanded.atoms`: Initial expanded LAMMPS configuration ($N = 7695$, $\rho = 0.001218$).
- `test_avbmc_reversibility.py`: Independent Python verification script proving microscopic reversibility (detailed balance) to machine precision.
- `analyze_avbmc_benchmark.py`: Analysis script computing cluster size distributions, kinetics, acceptance rates, and generating 4-panel publication plots.
- `avbmc_benchmark_comparison.png`: 4-panel high-resolution (300 DPI) comparison figure.
- `avbmc_reversibility_benchmark.tex`: Comprehensive publication-ready LaTeX article.
- `avbmc_reversibility_benchmark.pdf`: Compiled 7-page PDF documentation report.
- `run_with_avbmc/`: Full 200-sweep benchmark run with AVBMC enabled.
- `run_without_avbmc/`: Full 200-sweep canonical benchmark run without AVBMC.
- `submit_all.sh`: Master SLURM submission script for partition `gpu_v4_short`.

---

## Key Algorithms

### 1. Geometric Surface Asymmetry for Cluster Border Detection
To bias associations and dissociations specifically to cluster surfaces rather than dense cluster cores, each particle $i$ in a DBSCAN-detected cluster is evaluated using the geometric dipole asymmetry parameter:
$$\eta_i = \frac{1}{z_i} \left\| \sum_{j \in \text{neigh}(i)} \hat{\mathbf{r}}_{ij} \right\|$$
- **Cluster interior:** Near-isotropic neighbor distribution leads to vector cancellation ($\eta_i \approx 0$).
- **Cluster border:** Neighbors are confined to a hemisphere ($\eta_i \ge \eta_{\text{thresh}} = 0.5$).
Particles with $\eta_i \ge 0.5$ define the surface border pool $\mathcal{B}$.

### 2. Microscopic Reversibility (Detailed Balance)
AVBMC in-moves sample a free monomer from the gas phase into a spherical binding volume $V_{\text{bind}} = \frac{4}{3}\pi (1.5 r_{\text{cut}})^3$ centered at a surface target atom, while $K-1$ dummy trials sample the bulk box $V_{\text{box}}$ to compute the Rosenbluth weight $W_{\text{bulk}}$.

Out-moves sample $K$ trial positions in $V_{\text{box}}$ to extract a border particle into the gas phase.

The forward and reverse proposal arguments satisfy the exact mathematical identity:
$$\text{arg}_{\text{fwd}} \cdot \text{arg}_{\text{rev}} \equiv 1.000000000000$$
which guarantees:
$$\frac{P_{\text{acc}}(X \to Y)}{P_{\text{acc}}(Y \to X)} = \frac{e^{-\beta E(Y)}}{e^{-\beta E(X)}} \frac{\alpha(Y \to X)}{\alpha(X \to Y)}$$
Numerical testing with `test_avbmc_reversibility.py` across 5,000 independent transitions confirms detailed balance with:
- **Mean numerical discrepancy:** $3.1714 \times 10^{-17}$
- **Max numerical discrepancy:** $4.4409 \times 10^{-16}$ (machine epsilon)

### 3. $O(1)$ Host Neighbor-Cell Energy Optimization
During the $K$ Rosenbluth trials, energies are evaluated locally using an $O(1)$ 27-cell neighborhood query on the host CPU. Because cell occupancy is low ($\sim 0.03$ atoms/cell), trials take $< 1\,\mu\text{s}$ and completely bypass the PCIe bus. Overlapping trials ($r < 0.5\sigma$) are rejected immediately ($w = 0$). Grid repopulation and GPU synchronization occur only once per AVBMC step if moves are accepted.

---

## Benchmark Results (200 Sweeps on NVIDIA RTX PRO 4500 Blackwell)

| Metric | WITH AVBMC | WITHOUT AVBMC |
| :--- | :---: | :---: |
| **Total Particles $N$** | 7,695 | 7,695 |
| **Number Density $\rho^*$** | 0.001218 | 0.001218 |
| **Box Edge $L$** | 184.8576 | 184.8576 |
| **Reduced Temperature $T^*$** | 0.1540 | 0.1540 |
| **Final Cluster Count $N_{\text{clust}}$** | **60** | 67 |
| **Maximum Cluster Size $S_{\max}$ (Final / Peak)** | **1,099** / **1,374** | 1,063 / 1,068 |
| **Total Clustered Particles** | **7,631** | 7,545 |
| **Clustered Fraction (%)** | **99.17%** | 98.05% |
| **AVBMC In-Move Acceptance $P_{\text{AV, in}}$ (Avg / Final)** | **54.07%** / **55.62%** | N/A |
| **Canonical Translation Acceptance $P_{\text{Trans}}$** | 41.12% | 41.60% |
| **Canonical Rotation Acceptance $P_{\text{Rot}}$** | 55.48% | 55.51% |
| **Total GPU Execution Time on Blackwell** | 94.75 s | 92.10 s |
| **GPU Time per MC Sweep** | **473.7 ms** | 460.5 ms |
| **Energy Drift $\Delta E_{\text{drift}}$** | **0.0000E+00** | **0.0000E+00** |

---

## Reproducing the Verification and Benchmarks

1. **Verify Detailed Balance:**
   ```bash
   python3 test_avbmc_reversibility.py
   ```
2. **Submit SLURM Benchmarks:**
   ```bash
   ./submit_all.sh
   ```
3. **Analyze Outputs and Generate Plots:**
   ```bash
   python3 analyze_avbmc_benchmark.py
   ```
4. **Compile LaTeX Report:**
   ```bash
   module load texlive/20260310-GCC-13.2.0
   pdflatex avbmc_reversibility_benchmark.tex
   ```
