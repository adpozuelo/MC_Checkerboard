# MC_Checkerboard Examples and Benchmarks Suite

This directory contains reference simulation setups for **MC_Checkerboard**, demonstrating GPU-accelerated Checkerboard Monte Carlo, advanced cluster and identity swap algorithms, tabulated potential evaluations, and coupled Hybrid Monte Carlo (HMC) molecular dynamics trajectories via LAMMPS.

---

## Quick Reference Summary

| Directory | Model (`model`) | Physical System | Key Features & Algorithms | Execution Command |
| :--- | :--- | :--- | :--- | :--- |
| [`HS/`](file:///home/e.lomba/MC_Checkerboard/examples/HS) | `HS` | Binary Hard-Sphere Mixture | Athermal overlap rejection, contact distance matrix $\sigma_{ij}$, GPU identity swaps, virial pressure estimation | `../../bin/MCCB-gpu datos.nml 0` |
| [`LJ/`](file:///home/e.lomba/MC_Checkerboard/examples/LJ) | `LJ` | Binary Lennard-Jones Fluid | Continuous truncated & shifted 12-6 LJ, $\sigma_{ij}$ and $\varepsilon_{ij}$ matrices, GPU identity swap moves, $NVT$/$NpT$ | `../../bin/MCCB-gpu datos.nml 0` |
| [`LJG/`](file:///home/e.lomba/MC_Checkerboard/examples/LJG) | `LJG` | 7-Patch Colloidal Particles | Directional patchy attractions, angular and torsional constraints, quaternion rotations | `../../bin/MCCB-gpu input.d 0` |
| [`SSP/`](file:///home/e.lomba/MC_Checkerboard/examples/SSP) | `SSP` | Tetrahedral Patchy Particles | Palaia site-site patchy model, on-the-fly DBSCAN clustering, AVBMC aggregation moves | `../../bin/MCCB-gpu datos_ssp_tetrahedral.nml 0` |
| [`table_mixture/`](file:///home/e.lomba/MC_Checkerboard/examples/table_mixture) | `TABLE` | Binary Mixture ($N = 10,976$) | Tabulated Mie 50-49 interactions ($100,000$ points), GPU table evaluation, validated against LAMMPS energy | `../../bin/MCCB-gpu datos.nml 0` |
| [`HMC/`](file:///home/e.lomba/MC_Checkerboard/examples/HMC) | `SSP` | Tetrahedral Patchy Mixture | Hybrid MC: GPU Checkerboard MC coupled with LAMMPS rigid MD (`fix rigid/nve`), analytic & tabulated forces | `mpirun -np 1 ../../bin/MCCB-gpu datos_ssp_hmc.nml 0` |
| [`HMC_table_mixture/`](file:///home/e.lomba/MC_Checkerboard/examples/HMC_table_mixture) | `TABLE` | Binary Mixture ($N = 10,976$) | Hybrid MC: GPU Checkerboard MC coupled with LAMMPS GPU tabulated MD (`fix nve`, `pair_style table linear`) | `mpirun -np 1 ../../bin/MCCB-gpu datos.nml 0` |

---

## Detailed Case Descriptions

### 1. `HS/` — Non-Additive Hard-Sphere Binary Mixture
- **Physics**: Athermal hard-sphere mixture where configurations with core overlaps ($r_{ij} < \sigma_{ij}$) are strictly rejected ($P_{\text{acc}} = 0$), and non-overlapping configurations are unconditionally accepted ($P_{\text{acc}} = 1$).
- **Features**:
  - Symmetric contact distance matrix $\sigma_{ij}$.
  - GPU-accelerated checkerboard particle translations.
  - Non-destructive identity swap moves ($A \leftrightarrow B$) respecting hard-core non-additivity.
  - Virial equation of state pressure calculation from radial pair contact extrapolation.
- **Run**:
  ```bash
  cd examples/HS
  ../../bin/MCCB-gpu datos.nml 0
  ```

### 2. `LJ/` — Binary Lennard-Jones Fluid
- **Physics**: Classical continuous fluid mixture governed by 12-6 Lennard-Jones pairwise potentials:
  $$U_{ij}(r) = 4\varepsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - \left(\frac{\sigma_{ij}}{r}\right)^6\right] - U_{\text{shift}}$$
- **Features**:
  - Independent $\sigma_{ij}$ and $\varepsilon_{ij}$ interaction matrices.
  - Generalized GPU checkerboard identity swap moves with full Boltzmann energy difference checks.
  - Supports both canonical ($NVT$) and isobaric-isothermal ($NpT$) ensembles.
- **Run**:
  ```bash
  cd examples/LJ
  ../../bin/MCCB-gpu datos.nml 0
  ```

### 3. `LJG/` — Patchy Colloidal Particles (Lennard-Jones-Gauss)
- **Physics**: Anisotropic patchy colloids where spherical particles possess discrete attractive patches on their surface. Interactions combine isotropic Lennard-Jones repulsion/attraction with Gaussian angular modulation:
  $$V_{\text{patch}}(\mathbf{r}_{ij}, \mathbf{\Omega}_i, \mathbf{\Omega}_j) = V_{\text{radial}}(r_{ij})\, V_{\text{ang}}(\hat{\mathbf{r}}_{ij}, \hat{\mathbf{n}}_i, \hat{\mathbf{n}}_j)$$
- **Features**:
  - 7 patches per colloid with defined spatial unit vectors $\mathbf{p}_\alpha$ and symmetry axes $\mathbf{v}_{\text{op}}$.
  - Quaternions represent full 3D particle orientations.
  - Pairwise patch-patch coupling matrix $V_{\text{pot}}(\alpha, \beta)$ and angular widths $\sigma_{\text{ang}}$.
  - Trial rotations generated via random 3D rotation quaternions.
- **Run**:
  ```bash
  cd examples/LJG
  ../../bin/MCCB-gpu input.d 0
  ```

### 4. `SSP/` — Site-Site Patchy Colloids & Cluster Physics
- **Physics**: Tetrahedral patchy particle model developed by Palaia et al. Particles interact via an isotropic core repulsion with an attractive tail and 4 site-centered attractive spherical patches arranged in tetrahedral symmetry.
- **Features**:
  - Two-species binary patchy mixture ($N_{\text{species}} = 2$).
  - On-the-fly DBSCAN clustering engine identifying percolated networks, oligomers, and monomers.
  - Aggregation-Volume-Bias Monte Carlo (AVBMC) moves for efficient cluster condensation and evaporation across high free-energy barriers.
  - Geometric asymmetry analysis for cluster surface classification.
- **Run**:
  ```bash
  cd examples/SSP
  ../../bin/MCCB-gpu datos_ssp_tetrahedral.nml 0
  ```

### 5. `table_mixture/` — Tabulated Potential Binary Mixture
- **Physics**: Large-scale dense binary mixture ($N = 10,976$ particles) interacting via tabulated pair potentials derived from Mie 50-49 functional forms.
- **Features**:
  - $100,000$ points per table (`pot11.dat`, `pot12.dat`, `pot22.dat`), style `RSQ`.
  - Pair potentials are loaded into GPU texture/constant memory for high-throughput interpolation.
  - Initial configuration potential energy matches LAMMPS to $< 0.001\%$.
- **Run**:
  ```bash
  cd examples/table_mixture
  ../../bin/MCCB-gpu datos.nml 0
  ```

### 6. `HMC/` — Hybrid Monte Carlo for Patchy Colloids
- **Physics**: Couples spatial-decomposition Checkerboard Monte Carlo with microcanonical (NVE) molecular dynamics trajectory segments generated by LAMMPS.
- **Variants**:
  - `datos_ssp_hmc.nml`: Uses analytical rigid-body MD integration (`fix rigid/nve molecule`).
  - `datos_ssp_hmc_table.nml`: Uses tabulated force representations in both LAMMPS (`salr_lmp*.dat`) and MC GPU kernels (`salr_mc*.dat`).
- **Features**:
  - Seamless bidirectional memory transfer between Fortran/CUDA and LAMMPS C++ MPI library.
  - Rigorous Metropolis acceptance on total energy conservation $\Delta H = \Delta U + \Delta K$.
- **Run**:
  ```bash
  cd examples/HMC
  mpirun -np 1 ../../bin/MCCB-gpu datos_ssp_hmc.nml 0
  mpirun -np 1 ../../bin/MCCB-gpu datos_ssp_hmc_table.nml 0
  ```

### 7. `HMC_table_mixture/` — Hybrid Monte Carlo with Tabulated Mixtures
- **Physics**: Combines GPU Checkerboard MC with GPU-accelerated LAMMPS MD using tabulated Mie 50-49 potential tables for large systems ($N = 10,976$).
- **Features**:
  - Dynamic generation of LAMMPS table scripts (`pair_style table linear 100000`).
  - LAMMPS GPU acceleration package enabled (`-pk gpu 1 -sf gpu`).
  - Exact cell repopulation and GPU array re-synchronization upon accepted MD segments.
- **Run**:
  ```bash
  cd examples/HMC_table_mixture
  mpirun -np 1 ../../bin/MCCB-gpu datos.nml 0
  ```

---

## Namelist Configuration Reference

Input configuration files use Fortran 90 namelists divided into functional blocks:

### `&Control_Params`
- `istep_ini`: Starting sweep count (`0` for initial run, $> 0$ for restart).
- `istep_fin`: Total production sweeps.
- `Neq`: Equilibration sweeps discarded from statistics.
- `Nmove`: Trial move attempts per particle per sweep.
- `Nsave`: Frequency (in sweeps) for writing stdout and accumulating block averages.
- `Ndump`: Trajectory frame dump frequency.
- `Nrestart`: Checkpoint file save frequency (`data.restart`).
- `imovie`: Trajectory output toggle (`.true.` / `.false.`).
- `traj_format`: Format for trajectory (`'plain'`, `'netcdf'`, or `'both'`).
- `Npart_types`: Number of particle species in system.
- `lammps`: Master switch for Hybrid Monte Carlo LAMMPS MD moves (`.true.` / `.false.`, default: `.false.`). Only valid for `TABLE` and `SSP` models.
- `Ecpu_check`: Verifies initial and final system energy between CPU and GPU routines.

### `&MC_Params`
- `hmax`: Maximum translation displacement (in units of $\sigma$).
- `omax`: Maximum rotation displacement (in radians).
- `vmax`: Maximum box displacement components for $NpT$ volume changes.
- `displ_update`: Automatically adjusts `hmax`/`omax`/`vmax` to target $\sim 40\%$ acceptance.
- `temp0`, `temp1`: Initial and final reduced temperatures (supports linear thermal annealing).
- `npt`: Enables isobaric-isothermal ensemble ($NpT$).
- `pres`: Target reduced pressure $P^*$.
- `swap_moves`: Enables GPU checkerboard identity swap moves for multicomponent mixtures ($N_{\text{species}} \ge 2$) (`.true.` / `.false.`).
- `Nswap`: Number of swap sub-passes executed per swap attempt (default: 1).
- `Nswapf`: Cadence (in sweeps) of identity swap attempts (default: 1).

### `&Potential_Params`
- `model`: Potential interaction model: `'HS'`, `'LJ'`, `'LJG'`, `'SSP'`, or `'TABLE'`.
- `rangepp`: Radial potential cutoff in units of $\sigma$.
- `table_mc`: Enable tabulated potential evaluation in MC GPU kernels (`.true.` / `.false.`).
- `table_file_mc`: Prefix for MC potential table files (e.g., `'pot'` loads `pot11.dat`, `pot12.dat`, `pot22.dat`).

### `&Lammps_Params` *(Optional, active when `lammps = .true.`)*
- `timestep`: MD integration timestep in reduced time units (alias: `tstep_hmc`).
- `Nmd`: Number of MD integration steps per HMC trajectory segment (alias: `Nmd_hmc`).
- `hmc_freq`: Cadence (in sweeps) of HMC trial moves (alias: `Nmc_hmc`).
- `neigh_skin`: LAMMPS neighbor list skin thickness.
- `thermo_freq`: LAMMPS thermodynamic output frequency (`0` for quiet).
- `use_gpu`: Offloads LAMMPS pair styles to GPU (`-pk gpu 1 -sf gpu`).
- `table_lammps`: Force tabulated potential representation in LAMMPS.
- `table_file_lammps`: Prefix for LAMMPS potential table files.

---

## Environment Setup and Running Instructions

To execute simulations with OpenMPI and LAMMPS support on systems with NVIDIA HPC SDK:

```bash
# Load compiler and environment modules
source /home/e.lomba/bin/setup_trj

# Export OpenMPI and LAMMPS runtime library paths
export OPAL_PREFIX=/yggdrasil/usr_local/modules/x64_v4/software/NVHPC/25.3-CUDA-12.8.0/Linux_x86_64/25.3/comm_libs/12.8/hpcx/hpcx-2.22.1/ompi
export PATH=$OPAL_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$HOME/lammps_latest/build:$OPAL_PREFIX/lib:$LD_LIBRARY_PATH

# Run standard MC on GPU 0
../../bin/MCCB-gpu <input_file.nml> 0

# Run Hybrid MC with LAMMPS on GPU 0
mpirun -np 1 ../../bin/MCCB-gpu <input_file.nml> 0
```
