# GPU-Accelerated Monte Carlo for Patchy Particles

A high-performance GPU-accelerated Monte Carlo simulation code for patchy particle systems with anisotropic interactions in NVT and NpT ensembles.

## Authors

- **Eva González Noya** - eva.noya@iqf.csic.es
- **Antonio Díaz Pozuelo** - adiaz@iqf.csic.es

Instituto de Química Física Blas Cabrera (IQF-CSIC)

**Date:** February 2026

---

## Overview

This code implements a GPU-parallelized Monte Carlo simulation for systems of patchy particles with Lennard-Jones-Gauss (LJG) potential. The implementation uses a checkerboard cell decomposition scheme to enable conflict-free parallel Monte Carlo moves on GPU, achieving significant speedup compared to traditional CPU implementations.

### Key Features

- **GPU Acceleration**: CUDA Fortran implementation with checkerboard parallelization
- **Anisotropic Interactions**: Lennard-Jones core with angular (patch-patch) and torsional terms
- **Multiple Ensembles**: Supports both NVT (canonical) and NpT (isothermal-isobaric)
- **Flexible Particle Models**: Up to 7 patches per particle with customizable geometries
- **Adaptive MC Parameters**: Automatic adjustment of displacement parameters for optimal acceptance rates
- **Temperature Annealing**: Linear temperature ramping capability
- **Efficient Neighbor Lists**: Cell-based neighbor lists with periodic boundary conditions
- **Restart Capability**: Checkpoint and restart functionality for long simulations

---

## System Requirements

### Hardware
- NVIDIA GPU with CUDA Compute Capability 3.0 or higher
- Minimum 2 GB GPU memory (larger systems require more)

### Software
- NVIDIA HPC SDK (formerly PGI compiler suite) with CUDA Fortran support
- CUDA Toolkit 10.0 or higher
- Linux operating system (tested on CentOS/RHEL and Ubuntu)

---

## Code Structure

### Main Directory: `./src/`

#### Core Modules
- `Definitions.cuf` - Module definitions (precision, configuration, potential, properties, cell structure)
- `Main.cuf` - Main program with MC loop and simulation control

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
- `Read_input_data.cuf` - Read input parameters and configurations
- `Volume_move.cuf` - NPT volume change moves
- Various output routines for trajectories and properties

---

## Compilation

### Using the provided Makefile

```bash
mkdir bin
cd src
make
```

### Manual Compilation

```bash
nvfortran -O3 -Mcuda=cc70 -o ../bin/mc_gpu *.cuf
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

1. **`input.d`** - Simulation parameters
   - MC control (steps, moves, frequencies)
   - Temperature, pressure conditions
   - Potential model parameters
   - Particle type definitions
   - Patch geometries and interaction matrix

2. **`coords-[basename]`** - Initial configuration
   - Number of particles
   - Box matrix (h tensor)
   - Particle positions and types

### Running a Simulation

```bash
cd T0.154
../mc_gpu.exe
```

### Example Input Structure

See `T0.154/` directory for example input files.

### Output Files

- `movie.xyz` - Trajectory in XYZ format
- `movie-quat.dat` - Particle orientations (quaternions)
- `input-restart.d` - Restart input file
- `coords-[basename]-final` - Final configuration
- Standard output - Energy, acceptance rates, timing

---

## Simulation Parameters

### Monte Carlo Control
- `istep_ini`, `istep_fin` - Initial and final step numbers
- `Neq` - Equilibration steps
- `Nmove` - MC moves per sweep per particle
- `Nsave`, `Nsave2` - Output frequencies
- `Nrestart` - Restart file frequency

### Thermodynamic Parameters
- `temp0`, `temp1` - Initial and final temperatures (Kelvin)
- `pres` - Pressure (NpT ensemble)
- `npt` - Enable NpT ensemble (.true. / .false.)

### MC Displacement Parameters
- `hmax` - Maximum translation distance
- `omax` - Maximum rotation angle
- `vmax` - Maximum volume change (NpT)
- `displ_update` - Enable adaptive adjustment

### Potential Model
- **LJG (Lennard-Jones-Gauss)**: Supports multiple particle types with patches
  - `sigma_jon_aux` - Angular interaction width
  - `sigma_tor_jon` - Torsional interaction width
  - `rangepp` - Interaction cutoff

---

## Checkerboard Decomposition

The code uses a 3D checkerboard decomposition that divides the simulation box into cells:

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
2. **Configure input.d** with simulation parameters
3. **Equilibrate** at target temperature
4. **Production run** with trajectory output
5. **Analysis** of structural and thermodynamic properties

---

## Citation

If you use this code in your research, please cite:

*(Add publication reference here once available)*

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

---

## Acknowledgments

This work was supported by [funding agency/grant information].

Computing resources provided by [computing facility].

---

## Version History

- **v1.0** (February 2026) - Initial release
  - GPU-accelerated checkerboard Monte Carlo
  - LJG potential with patches
  - NVT and NpT ensembles
  - Comprehensive code documentation

