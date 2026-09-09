# Hybrid Monte Carlo (HMC) for Patchy Particle Systems (SSP Model)

This example demonstrates the Hybrid Monte Carlo (HMC) algorithm applied to tetrahedral patchy colloidal particles governed by the Site-Site Patchy (SSP) Palaia model. The simulation couples GPU checkerboard Monte Carlo with microcanonical molecular dynamics (MD) trajectory segments integrated via LAMMPS.

## System Description
- **Potential Model**: Site-Site Patchy (SSP) model by Palaia et al.
- **Topology**: Binary mixture of tetrahedral patchy particles ($N_{\text{species}} = 2$, 4 patches per particle arranged with tetrahedral symmetry).
- **Initial Configuration**: `data.beg` (LAMMPS format containing rigid body patch coordinates).
- **Ensemble**: Canonical ($NVT$) at $T^* = 0.15$.
- **Cluster Analysis**: On-the-fly DBSCAN cluster detection (`ncluster = 100`, $r_{\text{cl}} = 2.5$, $\text{minPts} = 7$).

## Available Configurations
1. **Analytic Potential Formulation (`datos_ssp_hmc.nml`)**:
   - LAMMPS uses rigid-body integration (`fix rigid/nve molecule`) with analytic potential formulations.
   - MC GPU routines evaluate interactions directly from analytical expressions.
2. **Tabulated Potential Formulation (`datos_ssp_hmc_table.nml`)**:
   - Tabulated force and energy profiles are used by both LAMMPS (`salr_lmp*.dat`) and the MC GPU kernel (`salr_mc*.dat`).
   - Uses `pair_style table linear` in LAMMPS.

## Execution
Run with MPI using OpenMPI and the compiled GPU executable:
```bash
# Environment setup
source /home/e.lomba/bin/setup_trj
export OPAL_PREFIX=/yggdrasil/usr_local/modules/x64_v4/software/NVHPC/25.3-CUDA-12.8.0/Linux_x86_64/25.3/comm_libs/12.8/hpcx/hpcx-2.22.1/ompi
export PATH=$OPAL_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$HOME/lammps_latest/build:$OPAL_PREFIX/lib:$LD_LIBRARY_PATH

# 1. Analytic SSP HMC
mpirun -np 1 ../../bin/MCCB-gpu datos_ssp_hmc.nml 0

# 2. Tabulated SSP HMC
mpirun -np 1 ../../bin/MCCB-gpu datos_ssp_hmc_table.nml 0
```
