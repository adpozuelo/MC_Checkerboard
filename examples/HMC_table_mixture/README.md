# Hybrid Monte Carlo (HMC) with Tabulated Potentials for Binary Mixtures

This example demonstrates the Hybrid Monte Carlo (HMC) implementation combining GPU-accelerated spatial-decomposition Checkerboard Monte Carlo with microcanonical (NVE) molecular dynamics trajectory segments executed in LAMMPS on the GPU, using user-supplied tabulated pair potentials.

## System Description
- **Components**: Binary fluid mixture ($N_{\text{species}} = 2$).
- **Total Particles**: $N = 10,976$ particles ($7,318$ of species 1 and $3,658$ of species 2).
- **Simulation Box**: Cubic box of edge length $L = 23.9396629$ (reduced density $\rho \approx 0.80$).
- **Interactions**: Tabulated pair potentials generated from Mie 50-49 functional forms:
  - $N_{\text{pts}} = 100,000$ points per table.
  - Potential files: `pot11.dat`, `pot12.dat`, `pot22.dat` (symlinked from `../table_mixture/`).
  - Cutoff radius: $r_{\text{cut}} = 1.20$.
- **Ensemble**: Canonical ensemble ($NVT$) at reduced temperature $T^* = 1.0$.

## Hybrid Monte Carlo Method
During the simulation:
1. Standard GPU checkerboard local translation moves are attempted in parallel across sub-lattices.
2. Every `hmc_freq` sweeps, an HMC trial move is initiated:
   - Atomic coordinates are transferred from CUDA device memory to LAMMPS.
   - Initial particle velocities are drawn from the Maxwell-Boltzmann distribution at $T^*$.
   - LAMMPS performs $N_{\text{md}}$ microcanonical (NVE) integration steps using `pair_style table linear 100000` accelerated via the LAMMPS GPU package (`-pk gpu 1 -sf gpu`).
   - Final coordinates and momenta are transferred back to the MC code.
   - The trial configuration is accepted or rejected according to the Metropolis-Hasting criterion on total energy conservation:
     $$\mathcal{P}_{\text{acc}} = \min\left(1,\, \exp\left[-\beta (\Delta U + \Delta K)\right]\right)$$
   - If accepted, the updated coordinates replace the current state; if rejected, the pre-MD state is restored.

## Configuration Parameters (`datos.nml`)
- **`&Control_Params`**:
  - `lammps = .true.`: Master flag enabling the LAMMPS HMC integration engine.
  - `Npart_types = 2`: Specifies a binary mixture.
  - `Ecpu_check = .true.`: Verifies energy consistency between CPU and GPU routines at start and termination.
- **`&Lammps_Params`**:
  - `timestep = 0.001`: MD integration timestep in reduced units ($\tau$).
  - `Nmd = 10`: Number of MD integration steps per HMC trajectory segment.
  - `hmc_freq = 5`: Frequency (in MC sweeps) of HMC moves.
  - `use_gpu = .true.`: Offloads LAMMPS pair force and neighbor list calculations to the GPU.
  - `table_lammps = .true.`: Instructs LAMMPS to use `pair_style table linear 100000`.
  - `table_file_lammps = 'pot'`: Table file prefix (`pot11.dat`, `pot12.dat`, `pot22.dat`).
- **`&Potential_Params`**:
  - `model = 'TABLE'`: Selects tabulated isotropic pair potential for MC GPU evaluations.
  - `table_mc = .true.`: Enables GPU table evaluation for checkerboard MC sweeps.
  - `table_file_mc = 'pot'`: Specifies table prefix for MC routines.

## How to Run
Ensure OpenMPI and LAMMPS shared libraries are in your environment:
```bash
# Load compiler and OpenMPI modules
source /home/e.lomba/bin/setup_trj
export OPAL_PREFIX=/yggdrasil/usr_local/modules/x64_v4/software/NVHPC/25.3-CUDA-12.8.0/Linux_x86_64/25.3/comm_libs/12.8/hpcx/hpcx-2.22.1/ompi
export PATH=$OPAL_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$HOME/lammps_latest/build:$OPAL_PREFIX/lib:$LD_LIBRARY_PATH

# Run HMC simulation on GPU device 0
mpirun -np 1 ../../bin/mc_gpu.exe datos.nml 0
```
