# ==============================================================================
# Makefile for MC_Checkerboard
# ==============================================================================
# The following environment variables can be defined or loaded via environment
# modules (e.g. 'module load NVHPC/25.3-CUDA-12.8.0 netCDF-Fortran/4.6.1-NVHPC-25.3-CUDA-12.8.0'):
#
# NVBIN       : path of nvfortran and MPI binaries
# NVINCLUDE   : path to nvfortran header files and modules
# NVLIBS      : path to nvfortran dynamic libraries
# NETCDFINC   : path to NETCDF and NETCDFF includes
# NETCDFLIB   : path to NETCDF and NETCDFF libraries
# LAMMPS_DIR  : root path to LAMMPS installation (defaults to ~/lammps_latest)
# LAMMPSINC   : path to LAMMPS fortran includes / mod files
# LAMMPSLIB   : path to LAMMPS dynamic libs (liblammps.so)
# FFTWINC     : path to FFTW3 include files (optional)
# FFTWLIB     : path to FFTW dynamic libs (optional)
# ==============================================================================

# LAMMPS path configuration (points to ~/lammps_latest)
LAMMPS_DIR    ?= $(HOME)/lammps_latest
LAMMPSINC      = $(LAMMPS_DIR)/fortran
LAMMPSLIB      = $(LAMMPS_DIR)/build

# Environment / Module fallbacks for NVHPC & NetCDF (cluster & local workstation)
NVHPC_ROOT     ?= $(wildcard /opt/nvidia/hpc_sdk/Linux_x86_64/2025 /opt/nvidia/hpc_sdk/Linux_x86_64/25.9 /opt/nvidia/hpc_sdk/Linux_x86_64/current)
NVHPC_DEFAULT   = $(firstword $(NVHPC_ROOT))

NVBIN         ?= $(if $(EBROOTNVHPC),$(EBROOTNVHPC)/Linux_x86_64/$(EBVERSIONNVHPC)/compilers/bin,$(if $(NVHPC_DEFAULT),$(NVHPC_DEFAULT)/compilers/bin))
NVINCLUDE     ?= $(if $(EBROOTNVHPC),$(EBROOTNVHPC)/Linux_x86_64/$(EBVERSIONNVHPC)/compilers/include,$(if $(NVHPC_DEFAULT),$(NVHPC_DEFAULT)/compilers/include))
NVLIBS        ?= $(if $(EBROOTNVHPC),$(EBROOTNVHPC)/Linux_x86_64/$(EBVERSIONNVHPC)/compilers/lib,$(if $(NVHPC_DEFAULT),$(NVHPC_DEFAULT)/compilers/lib))
NETCDFINC     ?= $(if $(EBROOTNETCDFMINFORTRAN),$(EBROOTNETCDFMINFORTRAN)/include,$(wildcard /usr/local/netcdf-nv/include))
NETCDFLIB     ?= $(if $(EBROOTNETCDFMINFORTRAN),$(EBROOTNETCDFMINFORTRAN)/lib,$(wildcard /usr/local/netcdf-nv/lib))
NETCDF_C_LIB  ?= $(if $(EBROOTNETCDF),$(EBROOTNETCDF)/lib,$(wildcard /usr/local/netcdf-nv/lib))

ifneq ($(NVBIN),)
export PATH   := $(NVBIN):$(PATH)
endif

# Compiler & Options (default to nvfortran; override default f77)
ifeq ($(origin FC),default)
FC             = nvfortran
endif
FC            ?= nvfortran
FCINC          = $(if $(NETCDFINC),-I$(NETCDFINC)) \
                 $(if $(NVINCLUDE),-I$(NVINCLUDE)) \
                 $(if $(FFTWINC),-I$(FFTWINC)) \
                 -I$(LAMMPSINC)

FCOPTS         = -O3 -cudalib=curand $(FCINC)

# Linker options and libraries (includes runtime rpath)
comma         := ,
RPATH_FLAGS    = -Wl$(comma)-rpath$(comma)$(LAMMPSLIB) \
                 $(if $(NETCDFLIB),-Wl$(comma)-rpath$(comma)$(NETCDFLIB)) \
                 $(if $(NETCDF_C_LIB),-Wl$(comma)-rpath$(comma)$(NETCDF_C_LIB))

LKLIBS         = $(if $(NETCDFLIB),-L$(NETCDFLIB)) \
                 $(if $(NETCDF_C_LIB),-L$(NETCDF_C_LIB)) \
                 $(if $(NVLIBS),-L$(NVLIBS)) \
                 $(if $(FFTWLIB),-L$(FFTWLIB)) \
                 -L$(LAMMPSLIB) -llammps -lnetcdff

LKOPTS         = -cuda -cudalib=curand $(RPATH_FLAGS)

# Output executable and directory
BINDIR         = ../bin
EXE            = $(BINDIR)/mc_gpu.exe

# Object files
OBJ            = Definitions.o \
                 potential_functions.o \
                 Read_input_data_nml.o \
                 Tools.o \
                 Subsweep_Energy_CUDA.o \
                 InitializeCheckerboard.o \
                 MCsweep_Checkerboard.o \
                 Shift_Cells.o \
                 CB_system_energy_gpu.o \
                 energy.o \
                 Volume_move.o \
                 UpdateCheckerboard.o \
                 netcdf_trajectory.o \
                 mod_avbmc.o \
                 mod_clusters.o \
                 lammps_hmc_wrapper.o \
                 Main.o

# Suffix rules
.SUFFIXES: .cuf .f90 .o

%.o: %.cuf
	$(FC) -c $(FCOPTS) $< -o $@

%.o: %.f90
	$(FC) -c $(FCOPTS) $< -o $@

.PHONY: all clean

all: $(EXE)

$(BINDIR):
	mkdir -p $(BINDIR)

$(EXE): $(OBJ) $(LAMMPSINC)/lammps.o | $(BINDIR)
	$(FC) $(LKOPTS) $(EXE_FLAG) -o $@ $(OBJ) $(LAMMPSINC)/lammps.o $(LKLIBS)

# Compile LAMMPS Fortran wrapper if needed
$(LAMMPSINC)/lammps.o: $(LAMMPSINC)/lammps.f90
	$(FC) -c $(FCOPTS) $< -o $@

# Module dependencies
potential_functions.o    : Definitions.o
Subsweep_Energy_CUDA.o   : Definitions.o potential_functions.o
Tools.o                  : Definitions.o potential_functions.o
Read_input_data_nml.o    : Definitions.o Tools.o
InitializeCheckerboard.o : Definitions.o Subsweep_Energy_CUDA.o Tools.o
MCsweep_Checkerboard.o   : Definitions.o Subsweep_Energy_CUDA.o Tools.o
Shift_Cells.o            : Definitions.o Subsweep_Energy_CUDA.o
CB_system_energy_gpu.o   : Definitions.o Subsweep_Energy_CUDA.o
energy.o                 : Definitions.o potential_functions.o Tools.o
Volume_move.o            : Definitions.o Subsweep_Energy_CUDA.o Tools.o
UpdateCheckerboard.o     : Definitions.o Subsweep_Energy_CUDA.o
netcdf_trajectory.o      : Definitions.o
mod_avbmc.o              : Definitions.o
mod_clusters.o           : Definitions.o mod_avbmc.o Subsweep_Energy_CUDA.o
lammps_hmc_wrapper.o     : Definitions.o $(LAMMPSINC)/lammps.o
Main.o                   : Definitions.o potential_functions.o Read_input_data_nml.o \
                           Tools.o Subsweep_Energy_CUDA.o InitializeCheckerboard.o \
                           MCsweep_Checkerboard.o Shift_Cells.o CB_system_energy_gpu.o \
                           energy.o Volume_move.o UpdateCheckerboard.o netcdf_trajectory.o \
                           mod_avbmc.o mod_clusters.o lammps_hmc_wrapper.o

clean:
	rm -f *.o *.mod $(EXE)
