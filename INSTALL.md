# Installation Guide for MC_Checkerboard

This guide provides step-by-step instructions for installing and building **`MC_Checkerboard`**, a high-performance, GPU-accelerated Monte Carlo simulation suite written in CUDA Fortran.

---

## 1. System Requirements

### Hardware Requirements
- **Host Architecture**: 64-bit x86 Linux (`x86_64`).
- **NVIDIA GPU**: CUDA Compute Capability 3.5 or higher (Compute Capability 7.0+ recommended: Tesla V100, A100, H100, RTX 30xx/40xx, or newer).
- **GPU Memory**: Minimum 2 GB VRAM (dense systems with $N \ge 10^4$ particles or large neighbor cell decompositions may require 4+ GB).

### Software Requirements
- **Compiler Suite**: **NVIDIA HPC SDK (NVHPC)** version 21.x or newer (tested with NVHPC 25.3), providing:
  - `nvfortran` (CUDA Fortran compiler)
  - CUDA Toolkit (CUDA 11.x or 12.x)
  - `cuRAND` random number generator runtime library
  - MPI compiler wrapper (OpenMPI / HPC-X)
- **NetCDF-Fortran & NetCDF-C**: Version 4.5+ (required for compressed binary `.nc` trajectory output).
- **LAMMPS (Optional / Recommended)**: Built as a shared library (`liblammps.so`) with the Fortran 2003 wrapper (`lammps.f90`), required for **Hybrid Monte Carlo (HMC)** moves (`model = 'TABLE'` and `'SSP'`).
- **FFTW3 (Optional)**: Required if using supplemental structural analysis or patch orientation tools.

---

## 2. Environment Setup

### 2.1 Loading Modules on HPC Clusters

On clusters using Environment Modules (Lmod or Tcl modules, such as EasyBuild installations), load the compiler and library modules:

```bash
# Load NVIDIA HPC SDK (provides nvfortran, CUDA, and cuRAND)
module load NVHPC/25.3-CUDA-12.8.0

# Load NetCDF C and Fortran libraries built with NVHPC
module load netCDF/4.9.2-NVHPC-25.3-CUDA-12.8.0
module load netCDF-Fortran/4.6.1-NVHPC-25.3-CUDA-12.8.0

# (Optional) Load OpenMPI / HPC-X and FFTW if not automatically loaded
module load FFTW/3.3.10-GCC-14.3.0
```

### 2.2 Setting Environment Variables Manually

If not using environment modules, ensure `nvfortran` and library paths are exported in your `~/.bashrc` or session:

```bash
# NVIDIA HPC SDK root (adjust path according to your installation)
export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/current
export PATH=$NVHPC_ROOT/compilers/bin:$PATH
export LD_LIBRARY_PATH=$NVHPC_ROOT/compilers/lib:$LD_LIBRARY_PATH

# NetCDF installation paths
export NETCDFINC=/usr/local/netcdf-nv/include
export NETCDFLIB=/usr/local/netcdf-nv/lib
export LD_LIBRARY_PATH=$NETCDFLIB:$LD_LIBRARY_PATH
```

---

## 3. Building LAMMPS for Hybrid Monte Carlo (Optional)

If you plan to use **Hybrid Monte Carlo (HMC)** trial moves (`lammps = .true.`), you must have a local installation of LAMMPS built as a shared library with its Fortran interface:

```bash
# 1. Clone LAMMPS repository
git clone -b release https://github.com/lammps/lammps.git ~/lammps_latest
cd ~/lammps_latest

# 2. Configure and build shared library
mkdir -p build && cd build
cmake ../cmake \
  -D BUILD_SHARED_LIBS=yes \
  -D PKG_MOLECULE=yes \
  -D PKG_RIGID=yes \
  -D PKG_GPU=yes \
  -D GPU_API=cuda \
  -D CMAKE_INSTALL_PREFIX=$HOME/lammps_latest

make -j4

# 3. Build Fortran 2003 wrapper module
cd ~/lammps_latest/fortran
nvfortran -c -fPIC lammps.f90 -o lammps.o

# 4. Add LAMMPS shared library to your runtime loader path
export LD_LIBRARY_PATH=$HOME/lammps_latest/build:$LD_LIBRARY_PATH
```

> **Note:** The `Makefile` in `MC_Checkerboard/src` defaults to looking for LAMMPS at `$(HOME)/lammps_latest`. If your LAMMPS installation resides elsewhere, set the `LAMMPS_DIR` environment variable before compiling:
> ```bash
> export LAMMPS_DIR=/path/to/your/lammps
> ```

---

## 4. Compiling `MC_Checkerboard`

1. Clone the `MC_Checkerboard` repository:
   ```bash
   git clone https://github.com/adpozuelo/MC_Checkerboard.git
   cd MC_Checkerboard
   ```

2. Enter the source directory and compile:
   ```bash
   cd src
   make clean
   make -j4
   ```

3. Upon successful compilation, the executable is created in the `bin/` directory:
   ```
   ../bin/MCCB-gpu
   ```
   *(Note: A backward-compatible symlink `../bin/mc_gpu.exe -> MCCB-gpu` is also created).*

### Makefile Configuration Variables

The build system in `src/Makefile` automatically discovers standard module paths (`$EBROOTNVHPC`, `$EBROOTNETCDFMINFORTRAN`, etc.) and system directories. You can customize the build by specifying or exporting the following variables:

| Variable | Description | Default |
| :--- | :--- | :--- |
| `FC` | Fortran compiler | `nvfortran` |
| `LAMMPS_DIR` | Root directory of LAMMPS source/build | `$(HOME)/lammps_latest` |
| `LAMMPSINC` | Path to LAMMPS Fortran module (`lammps.o`) | `$(LAMMPS_DIR)/fortran` |
| `LAMMPSLIB` | Path to `liblammps.so` | `$(LAMMPS_DIR)/build` |
| `NETCDFINC` | NetCDF header/module include directory | Auto-detected from `$EBROOTNETCDFMINFORTRAN` |
| `NETCDFLIB` | NetCDF-Fortran dynamic library directory | Auto-detected from `$EBROOTNETCDFMINFORTRAN` |
| `NVBIN` | Path to `nvfortran` binary directory | Auto-detected from NVHPC root |
| `NVINCLUDE` | NVHPC compiler include directory | Auto-detected from NVHPC root |
| `NVLIBS` | NVHPC compiler libraries directory | Auto-detected from NVHPC root |

Example of manual compilation invocation:
```bash
make -j4 LAMMPS_DIR=/opt/lammps NETCDFINC=/usr/include NETCDFLIB=/usr/lib64
```

---

## 5. Verifying the Installation

To verify that the executable was built correctly and can access the GPU:

### 5.1 Interactive Verification (on a GPU node or workstation)

```bash
cd ../examples/HS
../../bin/MCCB-gpu datos.nml 0
```

Expected terminal output:
```text
  ======================================================================
     Checkerboard Monte Carlo Simulation (CUDA Fortran Acceleration)
  ======================================================================
  -> Simulation Parameters initialized successfully.
  -> Total Particles: N = ...
  -> GPU Device initialized.
  ...
```

### 5.2 HPC Cluster Verification via Slurm

On HPC systems where frontend/login nodes lack GPUs, submit a short test batch job:

```bash
#!/bin/bash
#SBATCH --job-name=mc_test
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output=mc_test.out

# Load required modules
module load NVHPC/25.3-CUDA-12.8.0
module load netCDF/4.9.2-NVHPC-25.3-CUDA-12.8.0
module load netCDF-Fortran/4.6.1-NVHPC-25.3-CUDA-12.8.0

# Set runtime library paths
export LD_LIBRARY_PATH=$HOME/lammps_latest/build:$LD_LIBRARY_PATH

# Run test case
cd /path/to/MC_Checkerboard/examples/HS
../../bin/MCCB-gpu datos.nml 0
```

Submit with:
```bash
sbatch submit_test.sh
```

---

## 6. Troubleshooting & Common Issues

### Issue 1: `error while loading shared libraries: liblammps.so`
**Cause**: The dynamic linker cannot find `liblammps.so` at runtime.  
**Solution**: Ensure the path to your LAMMPS build directory is present in `LD_LIBRARY_PATH`:
```bash
export LD_LIBRARY_PATH=$HOME/lammps_latest/build:$LD_LIBRARY_PATH
```
Alternatively, check that the `-Wl,-rpath,...` flag in `src/Makefile` points to the correct directory.

### Issue 2: `error while loading shared libraries: libnetcdff.so`
**Cause**: The NetCDF-Fortran shared library is not loaded or not in the runtime library search path.  
**Solution**:
```bash
module load netCDF-Fortran
# or
export LD_LIBRARY_PATH=/path/to/netcdf/lib:$LD_LIBRARY_PATH
```

### Issue 3: `NVIDIA-SMI has failed because it couldn't communicate with the NVIDIA driver`
**Cause**: Executing `MCCB-gpu` directly on an HPC login/head node that does not possess an active NVIDIA GPU.  
**Solution**: Request an interactive GPU node (`salloc --gres=gpu:1 ...`) or submit via Slurm (`sbatch`).

### Issue 4: `NVFORTRAN-F-0004-Unable to open MODULE file cudafor.mod`
**Cause**: Compiler is unable to find CUDA Fortran intrinsic modules, typically because an unsupported or non-NVIDIA compiler (e.g. `gfortran`) was invoked instead of `nvfortran`.  
**Solution**: Verify that `which nvfortran` points to the NVIDIA HPC SDK and that `FC = nvfortran` is set.

### Issue 5: Adjusting Target GPU Architecture
By default, `-cuda` generates code for current compatible architectures. To target a specific GPU generation for optimal code generation, add `-gpu=ccXX` to `FCOPTS` in `src/Makefile` (e.g. `-gpu=cc70` for Volta V100, `-gpu=cc80` for Ampere A100, `-gpu=cc89` for RTX 4090, or `-gpu=cc90` for Hopper H100).
