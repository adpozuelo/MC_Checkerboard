#!/bin/bash
# ==============================================================================
# Script to compile LAMMPS as a shared library and build the Fortran wrapper
# dynamically on ladon27 using the NVHPC/CUDA environment.
# Supports incremental builds (resuming from last failure point).
# ==============================================================================
set -e

# Default parameters
CLEAN_BUILD=false
TARGET_DIR=""

# Parse arguments
for arg in "$@"; do
    if [ "$arg" == "--clean" ]; then
        CLEAN_BUILD=true
    else
        TARGET_DIR="$arg"
    fi
done

# Target directory argument (defaults to $HOME/lammps_latest)
LAMMPS_DIR="${TARGET_DIR:-$HOME/lammps_latest}"
# Expand relative paths to absolute paths
LAMMPS_DIR=$(eval echo "$LAMMPS_DIR")

echo "=========================================================================="
echo "          LAMMPS ON-THE-FLY BUILD SCRIPT (DYNAMIC ENVIRONMENT)            "
echo "=========================================================================="
echo "Target LAMMPS Directory: $LAMMPS_DIR"

# 1. Source user's compiler and module setup environment
SETUP_TRJ="$HOME/bin/setup_trj"
if [ -f "$SETUP_TRJ" ]; then
    echo "Sourcing environment from $SETUP_TRJ..."
    source "$SETUP_TRJ"
else
    echo "Warning: setup_trj environment script not found at $SETUP_TRJ."
    echo "Proceeding with current shell environment."
fi

# 1b. Bypass Binutils 2.42 ld crash (SIGILL) by preferring stable system ld (v2.35)
echo "Preferring system linker /usr/bin/ld over module binutils ld to bypass NVHPC linker crash..."
export PATH=/usr/bin:$PATH
hash -r

# 2. Verify compilers are available
if ! command -v nvfortran &> /dev/null; then
    echo "Error: nvfortran compiler not found. Please load the NVHPC module."
    exit 1
fi
if ! command -v mpicc &> /dev/null || ! command -v mpicxx &> /dev/null; then
    echo "Error: mpicc/mpicxx MPI compilers not found. Please load the MPI environment."
    exit 1
fi

# Dynamically locate the host compiler for CUDA (nvc++)
NVCXX_PATH=$(which nvc++)
if [ -z "$NVCXX_PATH" ]; then
    echo "Error: Could not locate nvc++ compiler in PATH."
    exit 1
fi
echo "Host CUDA Compiler resolved to: $NVCXX_PATH"

# 3. Clean up make-based files in source directory (required by CMake)
SRC_DIR="$LAMMPS_DIR/src"
if [ -d "$SRC_DIR" ]; then
    echo "Purging auto-generated files from make-based build system..."
    make -C "$SRC_DIR" purge || true
else
    echo "Error: Source directory $SRC_DIR does not exist."
    exit 1
fi

# 4. Set up build directory (incremental by default)
BUILD_DIR="$LAMMPS_DIR/build"
if [ "$CLEAN_BUILD" = true ]; then
    echo "Clean build requested. Recreating build directory at $BUILD_DIR..."
    rm -rf "$BUILD_DIR"
    mkdir -p "$BUILD_DIR"
else
    if [ ! -d "$BUILD_DIR" ]; then
        echo "Creating build directory at $BUILD_DIR..."
        mkdir -p "$BUILD_DIR"
    else
        echo "Reusing existing build directory for incremental compilation..."
    fi
fi
cd "$BUILD_DIR"

# 5. Configure CMake
echo "Configuring CMake..."
cmake ../cmake \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCUDA_HOST_COMPILER="$NVCXX_PATH" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_MPI=ON \
  -DBUILD_OMP=ON \
  -DBUILD_SHARED_LIBS=ON \
  -DPKG_EXTRA-PAIR=ON \
  -DPKG_GPU=ON \
  -DGPU_API=cuda \
  -DGPU_ARCH=auto \
  -DPKG_MOLECULE=ON \
  -DPKG_RIGID=ON \
  -DPKG_NETCDF=ON \
  -DNETCDF_DIR=/usr/local/netcdf-nv

# 6. Compile LAMMPS shared library
echo "Compiling LAMMPS dynamic library..."
cmake --build . -j 8

# 7. Compile Fortran Wrapper Interface
FORTRAN_DIR="$LAMMPS_DIR/fortran"
if [ -d "$FORTRAN_DIR" ]; then
    echo "Compiling Fortran wrapper module..."
    cd "$FORTRAN_DIR"
    nvfortran -O3 -c lammps.f90 -o lammps.o
else
    echo "Warning: Fortran directory $FORTRAN_DIR not found. Skipping wrapper build."
fi

echo "=========================================================================="
echo "                        BUILD COMPLETED SUCCESSFULLY                      "
echo "=========================================================================="
echo "Compiled shared library:  $BUILD_DIR/liblammps.so"
if [ -f "$FORTRAN_DIR/lammps.o" ]; then
    echo "Compiled Fortran wrapper: $FORTRAN_DIR/lammps.o"
    echo "Compiled module file:     $FORTRAN_DIR/liblammps.mod"
fi
echo "=========================================================================="
