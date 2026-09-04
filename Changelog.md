# Changelog

All notable changes to the **MC_Checkerboard** simulation code will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.7.0] - 2026-09-04

### Changed & Removed
- **Removed Mandatory MPI Runtime Requirement**:
  - Removed obsolete `Use mpi`, `MPI_Init`, and `MPI_Finalize` calls from `src/Main.cuf`.
  - The executable `mc_gpu.exe` can now be invoked directly as a standalone binary (e.g. `./bin/mc_gpu.exe datos.nml 0`) without needing `mpirun -np 1` or OpenMPI runtime daemon initialization.
  - Retained full backward compatibility for setups that continue to invoke `mpirun -np 1`.
- **Makefile & Build Modernization**:
  - Configured default Fortran compiler in `src/Makefile` to `FC ?= nvfortran` (overriding GNU Make's internal `f77` default).
  - Added robust auto-detection for local workstation NVHPC SDK installations (`/opt/nvidia/hpc_sdk/Linux_x86_64/2025`, `25.9`, `current`) and NetCDF paths (`/usr/local/netcdf-nv`).
  - Exported `PATH` from Makefile so `nvfortran` is resolved automatically without requiring manual environment module loads or shell exports on standard systems.
- **Launcher Scripts & Documentation**:
  - Updated `bin/run_MCGPU` and `examples/table_mixture/README.md` to execute `mc_gpu.exe` directly without OpenMPI wrappers.

---

## [2.6.0] - 2026-09-04

### Added
- **Tabulated Potential Mixtures (`model = 'TABLE'`)**:
  - Implemented GPU-accelerated Monte Carlo simulation for simple $n$-component fluid mixtures interacting via LAMMPS potential tables (`pot_int = 0`, `table_mc = .true.`).
  - Added support for both `RSQ` ($r^2$-uniform spacing from LAMMPS `pair_write ... rsq ...`) and linear `R` grid styles in `eval_table_pot_dev` and `eval_table_pot_host` (`src/potential_functions.cuf`).
  - Added 4-point cubic Catmull-Rom spline interpolation with boundary clamping on device constant memory and host.
  - Implemented flexible table loading in `load_table_potentials` (`src/Read_input_data_nml.cuf`):
    - Automatically discovers separate pairwise files (`pot11.dat`, `pot12.dat`, `pot22.dat`, or `mie11.dat`, `mie12.dat`, `mie22.dat`) or unified multi-table files.
    - Matches section keywords case-insensitively (`mie11`, `POT11`, `1_1`, `11`).
    - Robustly parses tokenized `N` line with style detection (`RSQ` vs `R`) and spacing calculation.
    - Dynamically configures the interaction cutoff matrix `rangep(i, j) = table_rmax(i, j)` directly from table headers, ensuring consistent checkerboard cell decompositions.
  - Integrated table evaluation into GPU checkerboard MC kernels:
    - `subsweep_LJ` for internal and neighbor cell moves with $r_{\min}$ overlap checks (`src/Subsweep_Energy_CUDA.cuf`).
    - Full-system and per-atom energy verification kernels `enerGPU_LJ` and `cell_per_atom_enerGPU_LJ`.
    - Host verification routines `ener` and `ener_LJ` (`src/energy.cuf`).
  - Created example benchmark case in `examples/table_mixture/` matching LAMMPS Mie 50-49 binary mixture data and verified that initial GPU energy matches LAMMPS to $< 0.001\%$.

### Fixed
- **OpenACC Dynamic Symbol Collision (`libgomp: TODO`)**:
  - Provided a `bind(C, name="acc_register_library")` stub in `src/Tools.cuf` to prevent NVHPC runtime from invoking GCC's unfinished `acc_register_library` implementation in `libgomp.so.1` (loaded by dynamic LAMMPS library dependencies).

---

## [2.5.1] - 2026-08-28

### Added
- **Adaptive Identity Swap Throttling**:
  - Automatically monitors rolling identity swap acceptance probability ($P_{\text{swap}}$) in `Main.cuf`.
  - When $P_{\text{swap}} < 0.00005$ due to high-density jamming, `Nswapf` is automatically throttled from 1 to 100 sweeps to avoid wasteful GPU kernel calls, and automatically restored to 1 if $P_{\text{swap}}$ recovers.
  - Added real-time color-coded terminal notices when throttling activates or restores.
- **Decoupled Move Frequencies (`Nvolf` and `Nswapf`)**:
  - Added `Nvolf` (NpT volume move frequency) and `Nswapf` (identity swap frequency) parameters to `&MC_Params` in `Read_input_data_nml.cuf` for decoupled sampling.
- **High-Accuracy Virial Pressure Calculation (Method A + Method B)**:
  - **Fine-Grain Log-Linear Contact Fit (Method A)**: Upgraded near-contact pair histogram in `src/energy.cuf` from 5 coarse $0.01\sigma$ shells to **10 fine shells of width $0.0025\sigma$** ($1.0000\sigma \to 1.0250\sigma$). Extrapolates contact correlation $g(\sigma_{ab}^+)$ via Log-Linear fitting ($\ln g(r) = a + b(r - \sigma_{ab})$) to match exponential near-contact peak in dense hard-sphere systems.
  - **Block Time-Accumulation (Method B)**: Added `mod_hs_virial_accum` module (`src/Definitions.cuf`) and `Accumulate_HS_Virial_Histogram` (`src/energy.cuf`), which periodically collects pair distance histograms across each `Nsave` window. Computes window-averaged $\langle g_{ab}(\sigma^+)\rangle$ and average box volume $\langle V \rangle$, reducing statistical noise by $\sim 3.16\times$ in NpT and NVT ensembles.

### Changed & Formatted
- **Terminal Progress Table Formatting**:
  - Aligned all header titles with exact field widths across active simulation feature sets.
  - Increased `P_swap` precision in progress tables to 5 decimal places (`f10.5`).

---

## [2.5.0] - 2026-08-27

### Added
- **GPU-Accelerated Identity Swap Moves ($A \leftrightarrow B$) for Binary Mixtures**:
  - Implemented CUDA device kernel `subsweep_HS_swap` for intra-cell swaps across all active cells in the 8 checkerboard subsets, concurrently testing the central cell and 26 neighbor cells across GPU warp threads in $O(1)$ time.
  - Implemented CUDA device kernel `subsweep_HS_cross_swap` for long-range cross-cell swaps between paired disjoint cells across opposite halves of the simulation volume ($\ge L/2 > 2\sigma_{\text{max}}$) with exact Hastings detailed balance:
    $$\alpha = \min\left(1, \frac{n_{A,1}\, n_{B,2}}{(n_{B,1} + 1)(n_{A,2} + 1)}\right)$$
  - Added host wrappers `subsweep_swap_gpu` and `subsweep_cross_swap_gpu` in `Subsweep_Energy_CUDA.cuf`.
  - Added driver subroutine `MCswap_Checkerboard` in `MCsweep_Checkerboard.cuf`.
  - Added host-device species type synchronization in `Convert_CB_RU` updating host `itype(i)` from device `List_Cell`.
  - Added namelist parameters `swap_moves` (`.true.` / `.false.`, default: `.false.`) and `Nswap` (sub-passes per MC cycle, default: `1`) to `&MC_Params` in `Read_input_data_nml.cuf`.
  - Added real-time tracking of swap acceptance rate (`P_Swap`) in main MC progress tables and cumulative 64-bit production counters (`nswap_att_tot`, `nswap_acc_tot`).
- **Binary Mixture & Model Consistency Checks**:
  - Implemented fatal error validation during namelist initialization stopping execution immediately if `swap_moves = .true.` is requested with `Npart_types /= 2` or non-HS models (`model /= 'HS'`).
  - Added internal safety guard returns (`If (.not. swap_moves .or. Npart_types /= 2 .or. pot_int /= -1) Return`) across all swap move routines.
- **Virial Pressure & Thermodynamics Enhancements**:
  - Integrated contact virial pressure accumulation into thermodynamic averages (`P_virial_av` in `Tools.cuf`).
  - Embedded virial pressure into per-atom stress tensor components in NetCDF trajectory frames (`netcdf_write_frame` in `netcdf_trajectory.cuf`).
  - Increased buffer cell safety factor in Hard Sphere volume moves (`Volume_move.cuf`) to prevent grid rebuild overflows during large volume expansions.

### Fixed
- **Contact Overlap Trap in Hard-Sphere Subsweep Kernel**:
  - Fixed strict strict inequality evaluation at line 2388 and 2411 of `Subsweep_Energy_CUDA.cuf` by updating contact boundary checks (`RSQ < 0.249999_wp`) to eliminate self-displacement false rejections on exact contact configurations.
- Suppressed repetitive cell info prints during periodic $NpT$ volume updates in `UpdateCheckerboard.cuf`.

### Documentation
- Updated `README.md` with algorithmic formulation, checkerboard parallelization mechanics, input parameters, and methodological citations for identity swap moves.
- Updated `examples/HS/datos.nml` with commented identity swap configuration options.

---

## [2.4.0] - 2026-08-20

### Added
- **Hard-Sphere (HS) Potential Model (`pot_int = -1`)**:
  - Athermal hard-sphere fluid and multicomponent mixture potential with support for non-additive cross diameters ($\sigma_{ij} \ne \frac{1}{2}(\sigma_{ii} + \sigma_{jj})$).
  - Dedicated `--- HS SIGMA MATRIX ---` block reading in `Read_input_data_nml.cuf`.
  - CPU verification routines in `energy.cuf` (`ener_HS`).
- **Contact Virial Pressure Calculation (`HS_Virial_Pressure`)**:
  - Multicomponent Lebowitz-Percus contact pair-correlation extrapolation from near-contact histogram shells.
  - Periodic `[P_HS]` progress output after equilibration.
- **Aggregation Volume Bias Monte Carlo (AVBMC)**:
  - Implemented AVBMC cluster association (bulk $\rightarrow V_{\text{in}}$) and dissociation ($V_{\text{in}} \rightarrow$ bulk) pair moves (`mod_avbmc.cuf`).
  - Configurable Rosenbluth (CBMC) $K$-trial scheme for bulk volume probing.
  - $O(1)$ single-particle cell list updates (`update_single_particle_cell`).
- **Progress Table Streamlining**:
  - Distinct equilibration (`Neq`) and production (`istep_fin`) phases with automated counter resets.
  - Consolidated single-line progress table with dynamic columns for cluster analysis and AVBMC acceptance.

---

## [2.3.1] - 2026-08-10

### Added
- **Cluster Border Points Analysis**:
  - Geometric asymmetry criterion based on normalized net neighbor displacement vectors.
  - Configurable threshold `asym_threshold` in `Control_Params` namelist.
  - Dedicated output trajectory files for clusters (`mclast_clconf.lammpstrj`) and border particles (`mclast_brdconf.lammpstrj`).

---

## [2.3.0] - 2026-08-01

### Added
- **LAMMPS Integration & Build System**:
  - Externalized LAMMPS potential definitions into modular `.lmp` include files (`potential_ssp_analytic.lmp`, `potential_ssp_table.lmp`).
  - Added compilation helper utility `utils/build_lammps.sh`.
  - Updated branding to `HYBRID MC_CB-GPU/LAMMPS SIMULATION ENGINE`.

---

## [2.2.0] - 2026-07-20

### Added
- **Tabulated Potential Interpolation**:
  - Support for prefix-based tabulated potentials with cubic spline interpolation on CPU and GPU.
  - GPU neighbor list integration (`neigh yes`) for LAMMPS HMC integration.
  - Static compile-time table memory structures eliminating dope vector lookups in CUDA kernels.

---

## [2.1.0] - 2026-07-10

### Added
- **Hybrid Monte Carlo (HMC)**:
  - Integrated multi-step LAMMPS Molecular Dynamics (MD) rigid-body segments with GPU checkerboard MC moves.
  - GPU acceleration for LAMMPS MD segments via `-pk gpu`.
  - Multi-site molecular mapping for patchy particles.
  - Energy per site (`En_tot/Ns`) logging in progress tables.

---

## [1.5.0] - 2026-06-25

### Added
- **GPU-Accelerated DBSCAN Cluster Analysis**:
  - Parallel BFS and neighbor searching kernels (`mod_clusters.cuf`).
  - $O(N)$ GPU cell-list neighbor finder.
  - Time-series cluster evolution export (`clusevol_mc.dat`).

---

## [1.4.0] - 2026-06-15

### Added
- **Site-Site Patchy (SSP / Palaia) Potential Model**:
  - Weeks-Chandler-Andersen (WCA) core with attractive cosine-squared tail.
  - Configurable scaling parameters: `sigp_factor`, `Rcp_factor`, `Rc_factor`, `epsp_factor`, and `patch_radial_factor`.
  - LAMMPS cross-validation scripts and test examples.
  - Custom input configuration file path via `data_file`.

---

## [1.3.0] - 2026-06-05

### Added
- Dedicated CUDA Subsweep LJ kernel for optimized single-particle translations.
- Mixed-precision CUDA implementation (`dkind` accumulator with single-precision arithmetic).
- Cell-based per-atom energy calculation routines on GPU.

---

## [1.2.0] - 2026-05-25

### Added
- Fortran namelist input configuration (`Read_input_data_nml.cuf`).
- Isotropic Lennard-Jones (LJ) fluids and multicomponent LJ mixtures.
- Parallel per-atom energy calculation on GPU.
- CPU energy validation check flag (`Ecpu_check`).

---

## [1.1.0] - 2026-05-15

### Added
- NetCDF binary trajectory export format (`netcdf_trajectory.cuf`) with AMBER conventions, cell dimensions, and quaternion metadata.

---

## [0.9.0] - 2026-05-01

### Added
- Initial 3D checkerboard cellular decomposition for GPU-accelerated Monte Carlo simulations.
- Support for patchy particles with Lennard-Jones-Gauss (LJG) potentials and Kern-Frenkel angular/torsional modulations.
- Canonical ($NVT$) and Isobaric-Isothermal ($NpT$) ensembles.
