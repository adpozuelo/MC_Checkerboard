# Sanchez-Burgos 50:50 Scaffold + Surfactant Mixture with AVBMC & Tabulated Potentials

This directory contains a complete simulation case for the **50:50 binary mixture of 4-patch scaffolds and 3-patch surfactants**, combining **tabulated potentials** (`table_mc = .true.`) and **Aggregation Volume Bias Monte Carlo (AVBMC)** (`avbmc = .true.`) on GPU.

---

## 1. Physical Model Specifications

- **Composition:** 50% scaffold molecules, 50% surfactant (client) molecules ($N = 2000$, 1000 scaffolds, 1000 surfactants).
- **Density & Box:** $\rho^* = N \sigma^3 / V = 0.1360$, cubic box of edge $L = 24.499865\sigma$.
- **Thermodynamic State:** Reduced temperature $T^* = 0.09$, thermostat temperature $T = 1.0$ ($k_B T = 1.0$, $\beta = 1.0$).

### Molecular Geometry & Patch Topology
- **Scaffold (Species 1, Type 0 in MC / Type 1 in LAMMPS):**
  - Central core diameter: $\sigma = 1.0$
  - 4 tetrahedral patches located at radial distance $r = \sigma/2 = 0.5\sigma$ from the center.
  - Normalized patch directions: $(\pm 1, \pm 1, \pm 1)/\sqrt{3}$.
  - Tetrahedral angle: $\arccos(-1/3) \approx 109.4712^\circ$.
- **Surfactant (Species 2, Type 1 in MC / Type 2 in LAMMPS):**
  - Central core diameter: $\sigma = 1.0$
  - 3 coplanar patches in a regular triangle located at radial distance $r = \sigma/2 = 0.5\sigma$ from the center.
  - Normalized patch directions: $(\cos(2\pi k/3), \sin(2\pi k/3), 0)$ for $k = 0, 1, 2$.
  - Inter-patch angle: $120.0^\circ$.

### Interactions
1. **Core-Core Interaction (All Pairs): Pseudo-Hard Sphere (PHS, Mie 50-49)**
   $$U_{\text{PHS}}(r) = \begin{cases} C \epsilon_R \left[ \left(\frac{\sigma}{r}\right)^{50} - \left(\frac{\sigma}{r}\right)^{49} \right] + \epsilon_R, & r < \left(\frac{50}{49}\right)\sigma \\ 0, & r \ge \left(\frac{50}{49}\right)\sigma \end{cases}$$
   with $\epsilon_R = 2/3$, $C = 50 (50/49)^{49} \approx 136.0888$, cutoff $r_c = \frac{50}{49}\sigma \approx 1.020408\sigma$.

2. **Patch-Patch Interaction: Continuous Square Well (CSW)**
   $$U_{\text{CSW}}(d) = \begin{cases} -\frac{\epsilon_{\text{CSW}}}{2} \left[ 1 - \tanh\left(\frac{d - r_w}{\alpha}\right) \right], & d < r_c \\ 0, & d \ge r_c \end{cases}$$
   with well width $r_w = 0.12\sigma$, steepness $\alpha = 0.005\sigma$, depth $\epsilon_{\text{CSW}} = 1/T^* = 1/0.09 \approx 11.1111$, cutoff $r_c = 0.20\sigma$.

3. **Specificity Rule ("Only unlike sites interact"):**
   - **Scaffold patch - Surfactant patch (unlike):** $V_{\text{pot}} = 1.0$ (CSW active)
   - **Scaffold patch - Scaffold patch (like):** $V_{\text{pot}} = 0.0$ (ZERO / non-interacting)
   - **Surfactant patch - Surfactant patch (like):** $V_{\text{pot}} = 0.0$ (ZERO / non-interacting)
   - **Core - Patch:** Non-interacting ($V_{\text{pot}} = 0.0$)

---

## 2. Tabulated Potentials (`table_mc = .true.`)

All interactions are evaluated through 100,001-point `RSQ`-spaced table files directly matching LAMMPS `potentials.table`:
- `pot11.dat`: Scaffold core - Scaffold core (`PHS`, $r \in [0.7, 1.17]\sigma$)
- `pot13.dat`: Scaffold core - Surfactant core (`PHS`, $r \in [0.7, 1.17]\sigma$)
- `pot33.dat`: Surfactant core - Surfactant core (`PHS`, $r \in [0.7, 1.17]\sigma$)
- `pot24.dat`: Scaffold patch - Surfactant patch (`CSW`, $r \in [10^{-6}, 1.17]\sigma$)
- `pot22.dat`: Scaffold patch - Scaffold patch (`ZERO`, $r \in [10^{-6}, 1.17]\sigma$)
- `pot44.dat`: Surfactant patch - Surfactant patch (`ZERO`, $r \in [10^{-6}, 1.17]\sigma$)

Cubic spline interpolation is used on both CPU and GPU device memory (`eval_table_pot_dev` / `eval_table_pot_host`), ensuring bitwise consistency between host and device.

---

## 3. AVBMC (Association Volume Bias Monte Carlo)

- **Clustering Engine:** G-DBSCAN (`mod_clusters.cuf`) with core contact cutoff $r_{\text{cl}} = 1.20\sigma$ and $\text{minPts} = 3$ (matching the paper's local aggregation criterion).
- **Surface Border Detection:** Geometric dipole asymmetry parameter $\eta_i = \frac{1}{z_i} \|\sum_j \hat{\mathbf{r}}_{ij}\| \ge 0.5$, identifying surface molecules of the condensed phase.
- **Rosenbluth Scheme:** Paired in-moves and out-moves with $K = 10$ trials sampling the bulk box volume $V_{\text{box}}$ and binding volume $V_{\text{bind}} = \frac{4}{3}\pi (1.5 r_{\text{cut}})^3$.
- **Frequency:** Evaluated every `ncluster = 5` sweeps, capped at `avbmc_max_trials = 50` attempts per step.

---

## 4. Directory Contents

- `datos.nml`: Simulation control namelist, topology, patch vectors, and $V_{\text{pot}}$ matrix.
- `data.atoms`: Initial 2000-molecule configuration (LAMMPS ellipsoid format).
- `setup_sb_avbmc.py`: Reproducible generator script to extract tables and convert geometries.
- `run_sb_avbmc.slurm`: SLURM batch job script targeting `gpu_v4_short`.
- `output.out`: Log of execution including energy checks and AVBMC acceptance rates.
- `pot*.dat`: Tabulated interaction files.

---

## 5. Running the Simulation

To submit on SLURM:
```bash
sbatch run_sb_avbmc.slurm
```

To run interactively on a GPU allocation:
```bash
/home/e.lomba/MC_Checkerboard/bin/MCCB-gpu datos.nml
```
