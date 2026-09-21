# Sanchez-Burgos 50:50 Mixture (3x3x3 Supercell, 54,000 Particles)

This directory contains the setup and simulation deck for the **3x3x3 supercell** of the Sanchez-Burgos 50:50 scaffold+surfactant mixture ($N = 54,000$ patchy ellipsoids) with **AVBMC** and **tabulated potentials** (`table_mc = .true.`).

---

## 1. System Specifications

- **Number of Particles:** $N = 54,000$ (27,000 4-patch scaffolds + 27,000 3-patch surfactants).
- **Box Dimensions:** Cubic box $L = 73.499592\sigma$ ($V \approx 397,060\sigma^3$).
- **Reduced Density:** $\rho^* = N \sigma^3 / V = 0.136000$ (conserved from the original $N = 2,000$ system).
- **Source Configuration:** Replicated from `../sb_avbmc_50_50/data.restart` via `replicate_restart.py`.
- **Interactions:** Tabulated Mie 50-49 PHS core + CSW patch attraction (`pot*.dat`).

---

## 2. Simulation Protocol (`datos.nml`)

- **Equilibration:** `Neq = 50000` steps.
- **Production:** `istep_fin = 100000` steps.
- **Total Sweeps:** 150,000 Monte Carlo sweeps.
- **AVBMC Moves:** Active (`avbmc = .true.`, `avbmc_k_trials = 10`, `avbmc_max_trials = 50`).
- **Cluster Analysis:** G-DBSCAN (`ncluster = 100`, `rcl = 1.20`, `minPts = 3`).
- **Checkpoints & Trajectory:** `Nrestart = 500`, NetCDF snapshot frequency `Ndump = 100`.

---

## 3. Short Queue Execution Feasibility

- **Benchmark Measured on GPU (`NVIDIA RTX PRO 4500 Blackwell` on `gpu_v4_short`):**
  - **0.057 seconds / MC sweep** for 54,000 particles (17.5 sweeps / second).
  - Total 150,000 sweeps runtime: $\approx 8,550\text{ s} \approx \mathbf{2.38\text{ hours}}$!
- **Partition Limit (`gpu_v4_short`):**
  - Default walltime limit: **3 days (72 hours)**.
  - Typical short walltime cutoff: **4 hours 15 minutes** or **24 hours**.
  - The job easily finishes in **$\sim$ 2.4 - 2.5 hours**, well within the short queue limit.

---

## 4. Execution

To run or submit to SLURM:
```bash
sbatch run_sb_avbmc_3x3x3.slurm
```
