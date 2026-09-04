# GPU-Accelerated Monte Carlo for Patchy Particles/Simple fluid mixtures

A high-performance GPU-accelerated Monte Carlo simulation code for patchy particle systems  with anisotropic interactionsi and simple fluid mixtures in NVT and NpT ensembles.

## Authors

- **Eva González Noya** - eva.noya@iqf.csic.es
- **Antonio Díaz Pozuelo** - adiaz@iqf.csic.es
- **Enrique Lomba** - enrique.lomba@csic.es

Instituto de Química Física Blas Cabrera (IQF-CSIC)

**Date:** August 2026

---

## Overview

This code implements a GPU-parallelized Monte Carlo simulation for systems of patchy particles with Lennard-Jones-Gauss (LJG) potential, Kern-Frenkel (K-F) and Lennard-Jones mixtures, as well as (possibly non-additive) Hard-Sphere (HS) mixtures. The implementation uses a checkerboard cell decomposition scheme [1] to enable conflict-free parallel Monte Carlo moves on GPU, achieving significant speedup compared to traditional CPU implementations.

### Note

Lennard-Jones reduced units are used. Potentials are truncated and shifted at rangepp*sigma_ij (note that other codes use by default a single cutoff)

### Key Features

- **GPU Acceleration**: CUDA Fortran implementation with checkerboard parallelization [1]
- **Multiple Potentials**: Hard-sphere (HS, including non-additive mixtures), Lennard-Jones (LJ) mixtures, tabulated potentials (LAMMPS potential tables for isotropic n-component mixtures), and patchy (LJG / SSP) models
- **GPU-Accelerated Identity Swaps**: Massively parallel identity swap moves ($A \leftrightarrow B$) for binary hard-sphere mixtures leveraging checkerboard cellular decomposition and long-range sublattice swaps to accelerate compositional mixing [5]
- **Anisotropic Interactions**: Lennard-Jones core with angular (patch-patch) and torsional terms
- **Multiple Ensembles**: Supports both NVT (canonical) and NpT (isothermal-isobaric)
- **Aggregation Volume Bias MC (AVBMC)**: Advanced cluster swap moves (association & dissociation) [3] with Configurable Rosenbluth (CBMC) bulk probing scheme and $O(1)$ single-particle cell list updates
- **Flexible Particle Models**: Up to 7 patches per particle with customizable geometries
- **Adaptive MC Parameters**: Automatic adjustment of displacement parameters for optimal acceptance rates
- **Temperature Annealing**: Linear temperature ramping capability
- **Efficient Neighbor Lists**: Cell-based neighbor lists with periodic boundary conditions
- **Hybrid Monte Carlo (HMC)**: Integration of GPU MC moves with LAMMPS Molecular Dynamics (MD) rigid body segments [4]
- **Restart Capability**: Checkpoint and restart functionality for long simulations
- **NetCDF Trajectory Output**: Efficient binary format with integrated quaternions and metadata (~35× smaller than text files)

---

## System Requirements

### Hardware
- NVIDIA GPU with CUDA Compute Capability 3.0 or higher
- Minimum 2 GB GPU memory (larger systems require more)

### Software
- NVIDIA HPC SDK (formerly PGI compiler suite) with CUDA Fortran support
- CUDA Toolkit 10.0 or higher
- Linux operating system (tested on CentOS/RHEL and Ubuntu)
- NetCDF-Fortran library (optional, for netCDF trajectory output)

---

## Code Structure

### Main Directory: `./src/`

#### Core Modules
- `Definitions.cuf` - Module definitions (precision, configuration, potential, properties, cell structure)
- `Main.cuf` - Main program with MC loop and simulation control
- `mod_avbmc.cuf` - Aggregation Volume Bias MC (AVBMC) module, CBMC Rosenbluth scheme, and border/bulk particle transfers

#### Energy and Interactions
- `energy.cuf` - CPU energy calculation for verification
- `CB_system_energy_gpu.cuf` - GPU energy calculation wrapper
- `Subsweep_Energy_CUDA.cuf` - CUDA kernels for energy and MC moves

#### Checkerboard Implementation
- `InitializeCheckerboard.cuf` - Setup cell decomposition and neighbor lists
- `UpdateCheckerboard.cuf` - Rebuild cells after volume changes
- `MCsweep_Checkerboard.cuf` - Perform MC sweeps over checkerboard sets
- `Shift_Cells.cuf` - Random cell boundary shifts to avoid artifacts

#### I/O and Utilities
- `Read_input_data_nml.cuf` - Read input parameters and configurations
- `Volume_move.cuf` - NPT volume change moves
- `netcdf_trajectory.cuf` - NetCDF trajectory output module
- Various output routines for trajectories and properties

---

## Compilation

### Prerequisites for NetCDF Support

If you want NetCDF trajectory output, load the NetCDF-Fortran module:

```bash
module load netCDF-Fortran  # or path to your netCDF installation
```

### Using the provided Makefile (check this file for environment variables that must be defined)

```bash
mkdir -p bin
cd src
make
```

The Makefile automatically includes NetCDF libraries if the proper environment variables are set.


### Manual Compilation

```bash
nvfortran -O3 -Mcuda=cc70 -o ../bin/mc_gpu.exe *.cuf
```

Adjust `-Mcuda=cc70` to match your GPU architecture:
- Tesla V100: `cc70`
- Tesla P100: `cc60`
- RTX 3090: `cc86`
- A100: `cc80`

---

## Usage

### Input Files

The simulation requires two input files:

1. **`filename.nml`** - Simulation parameters (1st argument of mc_cpu.exe)
   - MC control (steps, moves, frequencies)
   - Temperature, pressure conditions
   - Potential model parameters
   - Particle type definitions
   - Patch geometries and interaction matrix

2. **`data.atoms`** - Initial configuration (Lammps compatible `atom_style ellipsoid`)
   - Number of particles
   - Box matrix (h tensor)
   - Particle positions and types
   - Particle shape and quaternions

### Trajectory Output Formats

Configure output format in `next_input.nml`:
- **`plain`** - Text format (`movie.xyz`, `trajectory.lammpstrj`)
- **`netcdf`** - Binary NetCDF format (`trajectory.nc`) - **Recommended** (35× smaller files)
- **`both`** - Generate both formats

```
Trajectory format (plain, netcdf, both)
netcdf
```

### Running a Simulation

```bash
# Run using the default GPU device (ID 0):
mc_gpu.exe filename.nml

# Run using a specific GPU device (e.g., GPU 1):
mc_gpu.exe filename.nml 1
```

### Example Input Structure

Each potential model has a self-contained example directory under `examples/`, with its own namelist input and `data.atoms` configuration:

- **`examples/HS/`** - Non-additive binary hard-sphere mixture (`datos.nml`)
- **`examples/LJ/`** - Plain (isotropic) Lennard-Jones mixture (`datos.nml`)
- **`examples/LJG/`** - Lennard-Jones-Gauss patchy system with angular/torsional patches (`input.d`)
- **`examples/SSP/`** - Site-Site Patchy (tetrahedral, 4-patch) system (`datos_ssp_tetrahedral.nml`); also includes `data.atoms_lammps` / `data.atoms_lammps_reduced` and, at `examples/in.ssp`, `examples/forcefield.lj`, `examples/log.lammps`, the companion LAMMPS input/output used to cross-validate the SSP potential against LAMMPS (see [Mapping SSP to LAMMPS](#mapping-ssp-to-lammps-4) below)

The namelist filename passed as the first command-line argument to `mc_gpu.exe` is arbitrary (`datos.nml`, `input.d`, etc. are just naming conventions used across these examples).

### Output Files

**Trajectory Files** (format depends on `traj_format` selection):
- `movie.xyz` - Trajectory in XYZ format (if `plain` or `both`)
- `trajectory.nc` - NetCDF trajectory with coordinates, quaternions, cell data, and metadata (if `netcdf` or `both`)
- `trajectory.lammpstrj` - LAMMPS trajectory with coordinates, quaternions, cell data (if `plain` or `both`)

**Simulation Data**:
- `run-data.dat` - Time series of properties (energy, pressure, box dimensions)
- `data.restart` - Restart file configuration
- `input-restart.nml` - Restart input file with current parameters
- `clusevol_mc.dat` - Time evolution of clusters showing step, total clusters, maximum cluster size, total number of clustered particles, and percentage of clustered particles.
- `mclast_conf.lammpstrj` - Final full configuration output formatted in LAMMPS trajectory format for VMD visualization.
- `mclast_clconf.lammpstrj` - Final cluster configuration output formatted in LAMMPS trajectory format for VMD visualization.
- `mclast_brdconf.lammpstrj` - Final cluster border configuration based on geometric asymmetry criterion formatted in LAMMPS trajectory format.

**Standard Output**: Clean single-line progress table containing step index, total and per-site potential energies, translation/rotation acceptance ratios, moves per particle, cell grid size, CPU/GPU timing, and conditionally integrated cluster metrics (`N_Clust`, `Max_Cl`, `%Clust`) and AVBMC acceptance ratios (`P_AV_in`, `P_AV_out`). For **Hard-Sphere (HS)** runs only, an additional `[P_HS]` line reporting the instantaneous virial pressure is printed at the same cadence, once the equilibration phase (`Neq` steps) has finished (see [HS Virial Pressure](#virial-pressure-hs-only) below).

---

## Simulation Parameters (in namelist)

### Particle Type Indexing Convention

- **LAMMPS Data Files (`data.atoms`, `data.restart`) & Trajectories**: Follow standard 1-based LAMMPS particle type indexing (`1` to `Npart_types`). When reading `data.atoms`, 1-based particle types are converted to internal 0-based types (`itype = ityp - 1`).
- **Main Codebase & Internal Arrays**: Internal memory (`itype`) uses 0-based indexing (`0` to `Npart_types - 1`), where LAMMPS Type 1 corresponds to internal Type 0.
- **Input Namelists (`cluster_types`)**: Supports both 0-based (`0` to `Npart_types - 1`) and 1-based LAMMPS indexing (`1` to `Npart_types`) for user convenience. For example, both `cluster_types = 0` and `cluster_types = 1` target the first species (minority center particles).

### Monte Carlo Control
- `istep_ini`, `istep_fin` - Initial and final step numbers
- `Neq` - Equilibration steps
- `Nmove` - Number of trial moves per particle executed in each MC sweep (step). A single call to `MCsweep()` processes 100% of particles across all 8 checkerboard cell sets; every particle is subjected to `Nmove` trial moves (translations for isotropic models, or `Nmove/2` translations + `Nmove/2` rotations for anisotropic patchy models).
- `Nsave`, `Ndump` - Output frequencies
- `Nrestart` - Restart file frequency
- `data_file` - Path to the initial configuration file in LAMMPS format (default: `data.atoms`)
- `ncluster` - Perform cluster analysis every `ncluster` steps (default: `0`, which disables clustering)
- `rcl` - Cutoff distance for clustering (default: `1.8`). For patchy particle models, `rcl` should be set to the patch-patch interaction cutoff distance (e.g. $R_{\text{cut}} \approx 0.337$) to correctly identify physically bonded networks.
- `minPts` - Minimum number of neighbors for a core node in DBSCAN (default: `3`)
- `asym_threshold` - Geometric asymmetry parameter threshold (normalized net neighbor displacement vector magnitude) for identifying surface/border points of clusters (default: `0.5`)
- `cluster_types` - Integer array of up to 10 particle species to include in cluster analysis. Default is `-1` (unspecified). Accepts either 0-based (`0`) or 1-based (`1`) indexing for species selection (e.g., both `cluster_types = 0` and `cluster_types = 1` select the first species).

### Aggregation Volume Bias MC (AVBMC) Control
- `avbmc` - Enable AVBMC cluster association/dissociation moves (`.true.` / `.false.`, default: `.false.`)
- `border_criterion` - Surface/border particle detection method: `'energy'` (energy threshold relative to bulk) or `'asymmetry'` (geometric displacement asymmetry) (default: `'energy'`)
- `energy_border_ratio` - Energy ratio relative to bulk binding energy for border particle identification when `border_criterion = 'energy'` (default: `0.9`)
- `avbmc_k_trials` - Number of trial positions ($K$) generated in the bulk for Configurable Rosenbluth (CBMC) bulk volume probing (default: `10`)
- `avbmc_max_trials` - Maximum number of pair attempts per AVBMC step to bound execution time for large systems (default: `50`)

### Identity Swap Moves (Binary Mixtures)
- `swap_moves` - Enable GPU-accelerated identity swap moves ($A \leftrightarrow B$) (`.true.` / `.false.`, default: `.false.`).
  > **Note / Consistency Requirement:** Identity swap moves are **exclusively implemented for binary hard-sphere mixtures** (`Npart_types = 2` and `model = 'HS'`). Setting `swap_moves = .true.` with $N_{\text{part\_types}} \ne 2$ or non-HS models will trigger an immediate fatal error and cleanly terminate the simulation during input initialization.
- `Nswap` - Number of swap sub-passes executed per swap step (default: `1`). Each sub-pass performs 8 checkerboard subset intra-cell sweeps and 4 long-range cross-cell sweeps.
- `Nswapf` - Frequency (in MC sweeps) of identity swap moves (default: `1`, i.e., swap moves attempted every sweep when `swap_moves = .true.`).
- **Adaptive Swap Throttling**: When active, the main simulation controller monitors the rolling swap acceptance rate ($P_{\text{swap}}$). If $P_{\text{swap}}$ drops below $0.00005$ due to high-density compositional jamming, `Nswapf` is automatically throttled to `100` sweeps (saving GPU compute time) and automatically restored to `1` if $P_{\text{swap}}$ recovers.

### Thermodynamic Parameters
- `temp0`, `temp1` - Initial and final temperatures (Kelvin)
- `pres` - Pressure (NpT ensemble)
- `npt` - Enable NpT ensemble (`.true.` / `.false.`, default: `.false.`)
- `Nvolf` - Frequency (in MC sweeps) of NpT volume change moves (default: `1`, i.e., volume move attempted every sweep when `npt = .true.`)

### MC Displacement Parameters
- `hmax` - Maximum translation distance
- `omax` - Maximum rotation angle
- `vmax` - Maximum volume change (NpT)
- `displ_update` - Enable adaptive adjustment

### Potential Model
- **HS (Hard Sphere)**: Isotropic, athermal potential for hard-sphere fluids and **mixtures (not necessarily additive)**. Only the pair diameters `sigma_ij` are read (a full `Npart_types x Npart_types` matrix); there is **no epsilon matrix and temperature is not required** as input. Two particles overlap (forbidden configuration) when their center-center distance is below `sigma_ij`, otherwise the energy is zero. The Metropolis test reduces to a pure overlap test (a translation is accepted iff it creates no overlap). Because the mixture may be non-additive, `sigma_ij` can differ from `(sigma_ii + sigma_jj)/2` and is taken verbatim from the matrix. The contact distance `sigma_ij` also sets the checkerboard cell size, so no `rangepp` is needed.

#### Mathematical Formulation of HS
For the center-center distance $r$ between particles of types $i$ and $j$ with pair diameter $\sigma_{ij}$:
$$V_{ij}(r) = \begin{cases} \infty & r < \sigma_{ij} \\ 0 & r \ge \sigma_{ij} \end{cases}$$
The potential is athermal, so the temperature never enters the Metropolis acceptance in the NVT ensemble; a trial move is accepted if and only if it produces no overlap. For **non-additive** mixtures the cross diameter is independent, $\sigma_{ij} \ne \tfrac{1}{2}(\sigma_{ii}+\sigma_{jj})$, and is read directly from the `--- HS SIGMA MATRIX ---`.

> **NpT note:** In the NpT ensemble the acceptance rule is $\exp[-\beta P\,\Delta V + N\ln(V'/V)]$. Since HS is athermal the code fixes $\beta = 1$, so the input `pres` is interpreted as the **reduced pressure** $P^* = P/k_BT$ (equivalently $\beta P$).

#### Virial Pressure & Time-Accumulated Contact Extrapolation (HS only)

Because hard-sphere forces are impulsive, there is no smooth pairwise virial to sum as for LJ. Instead the code uses the exact Lebowitz-Percus contact theorem for a multicomponent hard-sphere mixture:
$$P = \rho T + \frac{2\pi}{3} T \sum_{a,b} \rho_a \rho_b\, \sigma_{ab}^3\, g_{ab}(\sigma_{ab}^+)$$
where $\rho_a = N_a/V$ is the partial number density of species $a$ and $g_{ab}(\sigma_{ab}^+)$ is the contact value of the pair correlation function for species pair $(a,b)$.

To achieve maximum precision even at high packing densities near the equation of state divergence, the Virial calculation combines two advanced schemes:
1. **Fine Shell Resolution & Log-Linear Contact Fit**: Pair separations are histogrammed into 10 fine shells of width $\delta_{ab} = 0.0025\,\sigma_{ab}$ right at contact ($r \in [\sigma_{ab},\, 1.025\,\sigma_{ab})$). The contact value $g_{ab}(\sigma_{ab}^+)$ is obtained via a Log-Linear least-squares fit ($\ln g(r) = a + b(r - \sigma_{ab})$), matching the physical exponential decay of pair correlations near contact.
2. **Block Time-Accumulation (NpT / NVT Compatible)**: Pair distance counts and box volume are accumulated periodically across all MC cycles within each `Nsave` window (`Accumulate_HS_Virial_Histogram`). The reported $P_{\text{virial}}$ is computed from the window-averaged contact distribution $\langle g_{ab}(\sigma_{ab}^+) \rangle$ and average volume $\langle V \rangle$, reducing statistical counting noise by $\sim 3.16\times$.

#### GPU Checkerboard Identity Swaps (Binary HS Mixtures)

In dense multicomponent fluid mixtures, traditional single-particle translation moves frequently encounter severe sampling bottlenecks caused by local steric cages (compositional jamming), resulting in sluggish structural relaxation and slow convergence. Identity swap moves ($A \leftrightarrow B$) [5] overcome this barrier by exchanging particle species identities without displacing atomic center-of-mass positions, dramatically accelerating phase space exploration and thermodynamic equilibration.

**Checkerboard Parallelization Strategy:**
1. **Intra-Cell Identity Swaps (`subsweep_HS_swap`)**:
   - In each of the 8 cellular checkerboard subsets, all active cells are separated by $\ge 1$ buffer cell from each other, ensuring completely conflict-free parallel execution.
   - For every active cell possessing both species ($n_A \ge 1$ and $n_B \ge 1$), a candidate pair $(a \in A, b \in B)$ is chosen.
   - GPU warp threads concurrently test the central cell and all 26 neighboring cells for hard-core overlaps under the swapped state ($a \to B, b \to A$).
   - If zero overlaps occur, the move is accepted unconditionally ($\alpha = 1.0$) because the cell composition $(n_A, n_B)$ remains invariant.
2. **Long-Range Cross-Cell Identity Swaps (`subsweep_HS_cross_swap`)**:
   - Spatially separated cell pairs $(C_1, C_2)$ situated at $(j_x, j_y, j_z)$ and $(j_x + N_{lx}/2, j_y + N_{ly}/2, j_z + N_{lz}/2)$ are paired across opposite halves of the box ($\text{separation} \ge L/2 > 2\sigma_{\text{max}}$).
   - Because their 26-neighborhoods are completely disjoint, both local neighborhoods are checked in parallel.
   - Acceptance satisfies exact detailed balance via the Hastings factor:
     $$\alpha = \min\left(1, \frac{n_{A,1}\, n_{B,2}}{(n_{B,1} + 1)(n_{A,2} + 1)}\right)$$
3. **Consistency & Constraints**:
   - Exclusively enabled for binary mixtures: `Npart_types = 2` and `model = 'HS'`. Attempting to use `swap_moves = .true.` with $N_{\text{part\_types}} \ne 2$ or continuous potentials triggers an immediate fatal error during initialization.
   - Preserves exact species stoichiometry ($N_A, N_B = \text{const}$).
   - Achieves sustained throughputs of $\sim 400,000\text{--}500,000$ swap attempts/second on modern GPUs with near-zero computational overhead.

#### Mathematical Formulation of LJ
For distance $r$:
$$V(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 - \text{shift\_offset} \right]$$
Where $\text{shift\_offset} = (\sigma/r_{\text{cutoff}})^{12} - (\sigma/r_{\text{cutoff}})^6$ shifts the potential to exactly $0.0$ at the cutoff $r_{\text{cutoff}} = rangepp \cdot \sigma$.

#### Mathematical Formulation of Tabulated Potentials (`TABLE`)
For general multicomponent mixtures interacting via custom or non-analytic pairwise potentials (such as Mie, Yukawa, EXP6, DPD, or numerical effective potentials generated by coarse-graining or LAMMPS `pair_write`), the code natively supports LAMMPS tabulated potential files:
$$V_{ij}(r) = \begin{cases}
V_{ij}(r_{\min}) & r < r_{\min} \\
\mathcal{S}_3\big( \{r_k, V_{ij,k}\}_{k=1}^N \big)(r) & r_{\min} \le r < r_{\max} \\
0 & r \ge r_{\max}
\end{cases}$$
Where $\mathcal{S}_3$ denotes a 4-point cubic Catmull-Rom spline interpolation evaluated directly on the GPU in constant device memory:
$$V(r) = \left(-\frac{1}{2}y_0 + \frac{3}{2}y_1 - \frac{3}{2}y_2 + \frac{1}{2}y_3\right) \xi^3 + \left(y_0 - \frac{5}{2}y_1 + 2y_2 - \frac{1}{2}y_3\right) \xi^2 + \left(-\frac{1}{2}y_0 + \frac{1}{2}y_2\right) \xi + y_1$$
Here $\xi = \text{idx\_real} - \lfloor\text{idx\_real}\rfloor$ is the fractional coordinate within the grid interval.

- **Grid Spacings**:
  - **`RSQ` Style**: Grid points are uniformly spaced in $r^2$ from $r_{\min}^2$ to $r_{\max}^2$ (standard output of LAMMPS `pair_write ... rsq ...` for optimal computational performance). Lookup index is:
    $$\text{idx\_real} = (r^2 - r_{\min}^2) \cdot \frac{N - 1}{r_{\max}^2 - r_{\min}^2}$$
  - **`R` Style**: Linear radial grid spacing:
    $$\text{idx\_real} = (r - r_{\min}) \cdot \frac{N - 1}{r_{\max} - r_{\min}}$$
- **Configuration in `datos.nml`**:
  ```fortran
  &Potential_Params
    model         = 'TABLE',   ! Selects tabulated potential mode
    table_mc      = .true.,    ! Enables MC table potential evaluations
    table_file_mc = 'pot'      ! File prefix (pot11.dat, pot12.dat, pot22.dat) or combined file
  /
  ```
- **Flexible File Discovery**: The loader automatically recognizes separate pairwise files with prefix `table_file_mc` (e.g. `pot11.dat`, `pot12.dat`, `pot22.dat` or `mie11.dat`), case-insensitive keywords (e.g. `mie11`, `POT11`, `1_1`), as well as single unified table files containing all pair sections.
- **Dynamic Interaction Range**: The cutoff for each pair $(i, j)$ is dynamically assigned from the table header $r_{\max}$, automatically establishing the checkerboard cell decomposition and interaction matrices without requiring manual `rangepp` or `sigma_LJ` definitions.

#### Mathematical Formulation of LJG
For colloid distance $r$:
$$E_{ij}(r, \Omega_i, \Omega_j) = V_{\text{LJ}}(r) \cdot \left[ F^-(r) + V_{\text{ang,tor}}^{\max}(\Omega_i, \Omega_j) F^+(r) \right]$$
Where:
1. **Lennard-Jones Potential**:
   $$V_{\text{LJ}}(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 - \text{shift\_offset} \right]$$
2. **Radial Switch Functions**:
   $$F^-(r) = \frac{1 - \tanh[a(r - \sigma_{op})]}{2}$$
   $$F^+(r) = \frac{1 + \tanh[a(r - \sigma_{op})]}{2}$$
   Here, $\sigma_{op} = xop \cdot \sigma$ is the distance where angular patch modulation switches on (with sharpness $a = 1000$). For $r < \sigma_{op}$, $F^-(r) \approx 1$ and $F^+(r) \approx 0$ (isotropic repulsion); for $r > \sigma_{op}$, $F^-(r) \approx 0$ and $F^+(r) \approx 1$ (anisotropic attraction).
3. **Angular and Torsional Modulation**:
   $$V_{\text{ang,tor}}^{\max}(\Omega_i, \Omega_j) = \max_{\alpha, \beta} \left\{ V_{\alpha\beta} \cdot V_{\text{ang}}(\theta_1, \theta_2) \cdot V_{\text{tor}}(\phi) \right\}$$
   * **Angular Alignment**:
     $$V_{\text{ang}}(\theta_1, \theta_2) = \exp\left( -\frac{\theta_1^2 + \theta_2^2}{\sigma_{ij}^{\text{ang}}} \right)$$
     Where $\theta_1$ (or $\theta_2$) is the angle between patch unit vector $\vec{p}_{i,\alpha}$ (or $\vec{p}_{j,\beta}$) and the intercolloid unit vector $\hat{r}_{ij}$ (or $-\hat{r}_{ij}$), and $\sigma_{ij}^{\text{ang}} = 2 \sigma_{\alpha}^{\text{ang}} \sigma_{\beta}^{\text{ang}}$ is the combined width.
   * **Torsional Alignment** (if `bool_tor == 1`):
     $$V_{\text{tor}}(\phi) = \max_k \left\{ \exp\left( -\frac{(\phi - \phi_{\alpha, k}^{\text{ref}})^2}{\sigma_{\text{tor}}} \right) \right\}$$
     Where $\phi$ is the dihedral angle between the reference orientation vectors of patch $\alpha$ and patch $\beta$ projected onto the perpendicular plane, and $\phi_{\alpha, k}^{\text{ref}}$ are allowed torsional alignment angles.

#### Mathematical Formulation of SSP

1. **Colloid Core-Core interaction**:
   For colloid-colloid center-center distance $r$:
   * $r \le R_c = 2^{1/6}\sigma$ (WCA Core):
     $$V_{\text{core}}(r) = -\epsilon_{\text{tail}} + 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 + \frac{1}{4} \right]$$
   * $R_c < r \le R_c \cdot Rc\_factor$ (Attractive Cosine-Squared Tail):
     $$V_{\text{core}}(r) = -\epsilon_{\text{tail}} \cos^2\left( \frac{\pi (r - R_c)}{2 (Rc\_factor - 1) R_c} \right)$$
   * $r > R_c \cdot Rc\_factor$:
     $$V_{\text{core}}(r) = 0$$

2. **Patch-Patch interaction**:
   For distance $d_p$ between patch $\alpha$ on particle $i$ and patch $\beta$ on particle $j$:
   * $d_p \le r_{\text{patch}}$:
     $$V_{\text{patch}}(d_p) = -V_{\alpha\beta}$$
   * $r_{\text{patch}} < d_p \le R_{cp}$ (Attractive Cosine-Squared Tail):
     $$V_{\text{patch}}(d_p) = -V_{\alpha\beta} \cos^2\left( \frac{\pi (d_p - r_{\text{patch}})}{2 (R_{cp} - r_{\text{patch}})} \right)$$
   * $d_p > R_{cp}$:
     $$V_{\text{patch}}(d_p) = 0$$
   Where $r_{\text{patch}} = sigp\_factor \cdot R_c$ is the patch radius, $R_{cp} = Rcp\_factor \cdot R_c$ is the patch interaction cutoff, and $V_{\alpha\beta}$ is the interaction strength defined in the `Vpot_matrix`.

#### Aggregation Volume Bias MC (AVBMC) & CBMC Rosenbluth Scheme

AVBMC moves [3] sample particle swaps between non-associated bulk particles and cluster border particles to accelerate cluster condensation and dissociation kinetics. To overcome low acceptance in dense bulk phases, a Configurable Rosenbluth (CBMC) scheme probabilistically probes bulk volume using $K$ trial positions:

1. **Dissociation Move ($V_{\text{in}} \rightarrow V_{\text{out}}$)**:
   When attempting to move a bound border particle from cluster binding volume $V_{\text{in}}$ to bulk $V_{\text{out}}$:
   - Generate $K$ independent trial positions $m=1 \dots K$ randomly in bulk $V_{\text{out}}$.
   - Evaluate interaction energy $u_m$ for each trial position on the GPU and calculate the Rosenbluth weight:
     $$W_{\text{new}} = \sum_{m=1}^K \exp(-\beta u_m)$$
   - Acceptance probability:
     $$\text{acc}(n \rightarrow o) = \min\left(1, \; \frac{N_{\text{out}} + 1}{N_{\text{in}} \cdot N_{\text{neigh}}} \cdot \frac{V_{\text{in}}}{V_{\text{out}}} \cdot \frac{W_{\text{new}}}{K \cdot \exp(-\beta u_{\text{cluster\_state}})}\right)$$
   - If accepted, select one of the $K$ bulk trial positions $m$ with probability $P(m) = \exp(-\beta u_m) / W_{\text{new}}$.

2. **Association Move ($V_{\text{out}} \rightarrow V_{\text{in}}$)**:
   When attempting to move a non-associated bulk particle into the binding shell $V_{\text{in}}$ of a cluster border particle:
   - Calculate Rosenbluth weight $W_{\text{old}}$ of the original bulk state using $K-1$ additional dummy bulk trial positions.
   - Sample a new target position in $V_{\text{in}}$ and compute its cluster interaction energy $u_{\text{cluster\_state}}$.
   - Acceptance probability:
     $$\text{acc}(o \rightarrow n) = \min\left(1, \; \frac{N_{\text{in}} \cdot N_{\text{neigh}} + 1}{N_{\text{out}}} \cdot \frac{V_{\text{out}}}{V_{\text{in}}} \cdot \frac{K \cdot \exp(-\beta u_{\text{cluster\_state}})}{W_{\text{old}}}\right)$$

#### Mapping SSP to LAMMPS [4]

To simulate this model in LAMMPS [4], we represent each colloid-patch assembly as a rigid molecule (with 1 center site + 4 patch sites) and configure a hybrid pair style:
* **Rigid Bodies**: Use `atom_style molecular` and define the rigid bodies via `fix rigid/nvt/small molecule` (or `fix rigid/small molecule`).
* **Pair Potential**: Use `pair_style hybrid/overlay` to overlay `lj/cut` and `cosine/squared` interactions.
* **WCA Shift**: To reproduce the $+1.0\epsilon$ shift in the core $V_{\text{core}}(r)$ equation, you must add **`pair_modify shift yes`** in the LAMMPS input script. This shifts the truncated `lj/cut` core interaction so that the energy goes smoothly to $0$ at the cutoff $R_c = 2^{1/6}\sigma$.
* **External Potential Scripts**: To keep the potential definition independent of the simulation wrapper, potential parameters are dynamically written at startup to external files and loaded into the LAMMPS execution stream via `include` statements:
  - `potential_ssp_analytic.lmp` is generated and included for site-site analytic hybrid overlay potentials.
  - `potential_ssp_table.lmp` is generated and included for tabulated potentials.


#### Implementing New Interactions

The codebase has a decoupled potential architecture that allows developers to easily add new interactions (CPU and GPU) in a single place:

1. **Namelist Input**: Add any new potential parameters to the `Potential_Params` namelist in [Read_input_data_nml.cuf](MC_Checkerboard/src/Read_input_data_nml.cuf).
2. **Mathematical Definition**: Define the potential math inside [potential_functions.cuf](MC_Checkerboard/src/potential_functions.cuf):
   * Implement a new subroutine declared with `attributes(device, host)` so it can compile for both CPU and GPU execution.
   * To keep it pure and thread-safe, pass all coordinate/orientation arrays and parameter matrices directly as dummy arguments (rather than referencing module globals).
3. **Routing**:
   * **CPU Routing**: In [energy.cuf](MC_Checkerboard/src/energy.cuf), update the master router `ener` to delegate calls to your new potential function based on `pot_int`.
   * **GPU Routing**: In [Subsweep_Energy_CUDA.cuf](MC_Checkerboard/src/Subsweep_Energy_CUDA.cuf), delegate calls inside the helper subroutines (such as `calc_pair_energy_sitesite` or `accumulate_cell_energy_sitesite`) to your new GPU device function.

---

## Checkerboard Decomposition

The code uses a 3D checkerboard decomposition [1] that divides the simulation box into cells:

- **8 Checkerboard Sets**: In 3D, cells are grouped into 8 sets based on (x,y,z) parity
- **Conflict-Free Parallelization**: All cells in a set can be updated simultaneously
- **Cell Size**: Automatically determined as ~1.1 × interaction range
- **Neighbor Lists**: Each cell tracks 26 neighbors with precomputed displacement vectors

This scheme enables efficient GPU parallelization while maintaining detailed balance.

---

## Performance Tips

1. **Cell Grid**: System should have at least 4×4×4 cells for optimal performance
2. **GPU Memory**: For large systems (>10,000 particles), ensure sufficient GPU memory
3. **Acceptance Rates**: Optimal ~40% for translations/rotations (automatically adjusted if `displ_update=.true.`)
4. **Equilibration**: Allow sufficient equilibration (typically 10⁴-10⁵ sweeps) before production

---

## Typical Workflow

1. **Prepare initial configuration** with desired density and particle orientations
2. **Configure filename.nml** with simulation parameters
3. **Equilibrate** at target temperature
4. **Production run** with trajectory output
5. **Analysis** of structural and thermodynamic properties

---

## Citation & References

If you use this code in your research, please cite:

```
Eva González Noya, Enrique Lomba, and Antonio Díaz Pozuelo, "GPU-Accelerated Monte Carlo for Simple Fluid and Patchy Particle mixtures: A high-performance GPU-accelerated Monte Carlo simulation code for patchy particle systems with anisotropic interactions and simple mixtures in NVT and NpT ensembles.", CSIC-Madrid (2026)
```

### Methodological References

1. **Checkerboard GPU Decomposition**:
   - J. A. Anderson, E. Jankowski, T. L. Grubb, M. Engel, and S. C. Glotzer, *"Massively parallel Monte Carlo for many-particle simulations on GPUs"*, *Journal of Computational Physics*, **254**, 27–38 (2013). DOI: [10.1016/j.jcp.2013.07.023](https://doi.org/10.1016/j.jcp.2013.07.023)

2. **Site-Site Patchy (SSP) Potential**:
   - I. Palaia and A. Šarić, *"Controlling cluster size in 2D phase-separating binary mixtures with specific interactions"*, *The Journal of Chemical Physics*, **156**, 194902 (2022). DOI: [10.1063/5.0087769](https://doi.org/10.1063/5.0087769)

3. **Aggregation Volume Bias Monte Carlo (AVBMC) Moves**:
   - T. D. Loeffler, A. Sepehri, and B. Chen, *"Improved Monte Carlo Scheme for Efficient Particle Transfer in Heterogeneous Systems in the Grand Canonical Ensemble: Application to Vapor–Liquid Nucleation"*, *Journal of Chemical Theory and Computation*, **11**(2), 542–552 (2015). DOI: [10.1021/ct500958y](https://doi.org/10.1021/ct500958y)
   - B. Chen and J. I. Siepmann, *"Aggregation-Volume-Bias Monte Carlo for Simulating Association in Physical Systems"*, *The Journal of Physical Chemistry B*, **104**(36), 8725–8733 (2000). DOI: [10.1021/jp001952u](https://doi.org/10.1021/jp001952u)

4. **LAMMPS Molecular Dynamics Engine**:
   - A. P. Thompson, H. M. Aktulga, R. Berger, D. S. Bolintineanu, W. M. Brown, P. S. Crozier, P. J. in 't Veld, A. Kohlmeyer, S. G. Moore, T. D. Nguyen, R. Shan, M. J. Stevens, J. Tranchida, C. Trott, and S. J. Plimpton, *"LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales"*, *Computer Physics Communications*, **271**, 108171 (2022). DOI: [10.1016/j.cpc.2021.108171](https://doi.org/10.1016/j.cpc.2021.108171)
   - S. Plimpton, *"Fast Parallel Algorithms for Short-Range Molecular Dynamics"*, *Journal of Computational Physics*, **117**, 1–19 (1995). DOI: [10.1006/jcph.1995.1039](https://doi.org/10.1006/jcph.1995.1039)

5. **Identity Swap Monte Carlo Moves**:
   - D. A. Kofke and E. D. Glandt, *"Monte Carlo simulation of multicomponent equilibria in a semigrand canonical ensemble"*, *Molecular Physics*, **64**(6), 1105–1131 (1988). DOI: [10.1080/00268978800100761](https://doi.org/10.1080/00268978800100761)
   - A. J. Schultz and D. A. Kofke, *"Semigrand canonical Monte Carlo simulation: Comparison of methods for identity exchange"*, *The Journal of Chemical Physics*, **133**(10), 104101 (2010). DOI: [10.1063/1.3486085](https://doi.org/10.1063/1.3486085)
   - D. Frenkel and B. Smit, *"Understanding Molecular Simulation: From Algorithms to Applications"*, 2nd ed., Academic Press, San Diego (2002).

---

## License

GNU GENERAL PUBLIC LICENSE Version 3

---

## Troubleshooting

### Common Issues

**"Your system is too small, sorry!"**
- Increase box size or decrease interaction range
- Ensure at least 4 cells per dimension

**"You exceeded the maximum number of Cells"**
- Reduce box size or increase `NCell_max` in code and recompile

**Low acceptance rates**
- Enable `displ_update=.true.` for automatic adjustment
- Manually reduce `hmax`, `omax`, or `vmax`

**GPU memory errors**
- Reduce system size
- Use GPU with more memory

### Getting Help

For questions or issues, contact the authors:
- Eva González Noya: eva.noya@iqf.csic.es
- Antonio Díaz Pozuelo: adiaz@iqf.csic.es
- Enrique Lomba: enrique.lomba@csic.es

---

## Acknowledgments

This work was supported by CSIC.

Computing resources provided by CSIC.

---

## Version History

For a complete record of all versions and features, see [Changelog.md](Changelog.md).

- **V2.7** (September 2026) Standalone Direct Execution & Modernized Build
  - Removed obsolete MPI runtime dependency (`MPI_Init`/`MPI_Finalize`) from `src/Main.cuf`, allowing `mc_gpu.exe` to be invoked directly from the command line without `mpirun -np 1`.
  - Updated `Makefile` to use `FC = nvfortran` with automatic detection for standard NVHPC and NetCDF library paths.

- **V2.6** (September 2026) Tabulated Potential Fluid Mixtures (`model = 'TABLE'`)
  - Full GPU-accelerated Monte Carlo support for $n$-component mixtures interacting via LAMMPS potential tables (`RSQ` and linear `R` styles) with Catmull-Rom cubic spline interpolation.
  - Flexible file discovery and keyword loader for multi-component interaction tables.
  - Validated against LAMMPS Mie 50-49 binary fluid mixture benchmark with $< 0.001\%$ energy agreement.

- **V2.5** (August 2026) GPU-Accelerated Identity Swaps, Non-Additive Hard-Sphere Potential & Virial Pressure
  - **GPU Checkerboard Identity Swaps**: Massively parallel identity swap moves ($A \leftrightarrow B$) for binary mixtures using intra-cell checkerboard warp evaluation (`subsweep_HS_swap`) and long-range disjoint cross-cell swaps (`subsweep_HS_cross_swap`) with exact Hastings detailed balance.
  - **Sanity Checks & Consistency**: Strict input validation enforcing `Npart_types = 2` and `model = 'HS'` with fatal error termination on unsupported configurations.
  - **Contact Overlap Trap Fix**: Resolved exact-contact self-displacement rejection in Hard Sphere CUDA subsweep kernel.
  - **Hard-Sphere (HS) Potential & Virial Pressure**: Athermal HS potential (`pot_int = -1`) for additive/non-additive mixtures, multicomponent contact virial pressure extrapolation (`HS_Virial_Pressure`), and NetCDF stress embedding.
  - **Namelist Parameters**: Added `swap_moves` and `Nswap` to `&MC_Params` and updated `examples/HS/datos.nml`.

- **V2.4** (July/August 2026) Aggregation Volume Bias Monte Carlo (AVBMC) & CBMC Rosenbluth Scheme
  - Implemented AVBMC cluster association (bulk $\rightarrow V_{\text{in}}$) and dissociation ($V_{\text{in}} \rightarrow$ bulk) pair moves (`mod_avbmc.cuf`).
  - Integrated Configurable Rosenbluth (CBMC) $K$-trial scheme for probabilistic bulk volume sampling ($W_{\text{new}} = \sum_{m=1}^K \exp(-\beta u_m)$).
  - Added namelist control parameters `avbmc`, `border_criterion`, `energy_border_ratio`, `avbmc_k_trials` (default `10`), and `avbmc_max_trials` (default `50`).
  - Implemented $O(1)$ single-particle cell list updates (`update_single_particle_cell`) with orientation preservation and global grid fallback search.
  - Corrected 0-indexed position array `r` slicing offsets (`(id-1)*ndim : (id-1)*ndim+2`) and replaced boundary wrapping with exact `Floor` functions, achieving 53%+ AVBMC move acceptance.
  - Streamlined main Monte Carlo progress table into a unified single-line format with dynamic cluster and AVBMC metrics.

- **V2.3.1** (July 2026) Cluster Border Points Analysis & Geometric Asymmetry Criterion
  - Implemented geometric asymmetry criterion (normalized net neighbor displacement vector magnitude) to identify surface/border particles in clusters.
  - Added configurable threshold `asym_threshold` (default `0.5`) in `Control_Params` namelist.
  - Renamed final output configuration files to `mclast_conf.lammpstrj`, `mclast_clconf.lammpstrj`, and `mclast_brdconf.lammpstrj` to avoid overwriting outputs when running post-processing utilities like `trj_analysis`.

- **V2.3** (August 2026) Multi-state (S1-A1 to S1-A3, S1-A1 to S2-A2, S2-A2 to S3-A3) and Hybrid (MC to HMC) transitions
  - LAMMPS potential definitions are now stored in external files `potential_ssp_analytic.lmp` and `potential_ssp_table.lmp` and included in the LAMMPS input script.
  

- **V2.2** (July 2026) Table Potential Interpolation and HMC Performance Improvements
  - Implemented prefix-based table potential loading and linear interpolation on both host (CPU) and device (GPU constant/global memory lookup).
  - Modified LAMMPS HMC initialization to use GPU neighbor lists (`neigh yes`) when table potentials are active, achieving a 35%+ speedup in HMC integration.
  - Added CPU/GPU timer measurements and total MD step logging to the HMC acceptance periodic printouts.
  - Reset periodic performance timers after each HMC step to keep MC sweep statistics clean.
  - Refactored table memory structure from dynamic allocation to fixed compile-time size, completely removing dope vector memory lookups inside the GPU lookup kernel.
  - Decoupled LAMMPS force field setups from the Fortran wrapper into separate `.lmp` potential script files (`potential_ssp_table.lmp` and `potential_ssp_analytic.lmp`) loaded dynamically via `include` statements.

- **V2.1** (July 2026) Hybrid Monte Carlo (HMC) strategy and Energy per Site outputs
  - Integrated HMC combining checkerboard GPU MC moves with multi-step LAMMPS Molecular Dynamics (MD) rigid body segments.
  - Enabled GPU targeting for HMC LAMMPS MD segments (via program command line argument) using the LAMMPS GPU package (`-pk gpu`).
  - Implemented exact molecular site-to-particle mapping for HMC SSP systems, treating each colloid particle (center + patches) as a single multi-site LAMMPS molecule.
  - Configured core WCA potential energy shifting (`pair_modify shift yes`) in LAMMPS for exact energy alignment with the MC engine.
  - Added the `En_tot/Ns` (energy per site) column to the periodic terminal output for direct comparison against LAMMPS `E_pair`.
  - Fixed out-of-bounds `itype` array indexing in HMC wrappers.

- **V1.6** (July 2026) Documentation updates for cluster analysis default configurations
  - Clarified targeting of interaction patches (external types 2 and 4, internally 1 and 3) during DBSCAN cluster search.
  - Documented physical setup of the `rcl` cutoff distance using the patch-patch interaction range.

- **V1.5** (July 2026) GPU-accelerated DBSCAN Cluster Analysis
  - Added GPU-accelerated DBSCAN cluster identification using neighbor search and parallel BFS kernels.
  - Periodic cluster analysis triggered via namelist configuration (`ncluster`, `rcl`, `minPts`, `cluster_types`).
  - Output metrics printed to stdout and saved to `clusevol_mc.dat`.
  - Resolved background cell-grid energy accumulation bug.

- **V1.4** (June 2026) Configurable patchy potential scaling and bug fixes
  - Added namelist potential scaling parameters (`sigp_factor`, `Rcp_factor`, `Rc_factor`) for the Palaia site-site potential model.
  - Relocated initialization sanity checks to prevent Segmentation Fault during startup.

- **V1.3** (June 2026) Support for LJ mixture
  - Removed unneeded moves for simple systems
  - Output info adapted to LJ mixture

- **v1.2** (May 2026) Added compatibility with trj_analysis tool.
  - Input files transformed to namelist
  - Configuration files converted to LAMMPS system.data style
  - Output information beautified

- **v1.1** (February 2026) - NetCDF trajectory support
  - Added `netcdf_trajectory.cuf` module for binary trajectory I/O
  - Configurable trajectory format via `traj_format` parameter
  - AMBER-compatible NetCDF structure with quaternions
  - 35× file size reduction vs plain text
  - Verification utilities included
  
- **v1.0** (February 2026) - Initial release
  - GPU-accelerated checkerboard Monte Carlo
  - LJG potential with patches
  - NVT and NpT ensembles
  - Comprehensive code documentation

