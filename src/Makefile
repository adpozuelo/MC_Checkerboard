# ==============================================================================
# Makefile for MC_Checkerboard
# ==============================================================================
# Targets:
#   make          build both host-architecture binaries + launcher (default)
#   make v3       build only bin/MCCB-gpu-v3 (-tp=x86-64-v3: AVX2, e.g. ladon25/26)
#   make v4       build only bin/MCCB-gpu-v4 (-tp=x86-64-v4: AVX-512)
#   make clean    remove build/, binaries, launcher and legacy objects
# bin/MCCB-gpu is a launcher that runs the v4 binary on CPUs with avx512f and
# the v3 binary otherwise (bin/mc_gpu.exe is a compatibility symlink to it).
# Objects/modules go to build/<arch>; GPU code targets cc75,86,89,90,120.
#
# Overridable variables (environment or command line):
#   NVBIN / NVINCLUDE / NVLIBS : nvfortran bin, include and lib directories
#   NETCDFINC / NETCDFLIB / NETCDF_C_LIB / HDF5LIB : NetCDF(-Fortran)/HDF5 paths
#                                (default on Ladon: /usr/local/modules/x64_<arch>/...)
#   LAMMPS_DIR (~/lammps_latest) : LAMMPS root (LAMMPSINC/LAMMPSLIB derive from it;
#                                env values are ignored, override on the command line)
#   FFTWINC / FFTWLIB          : FFTW3 include/lib (optional)
# Library paths are embedded as DT_RPATH (--disable-new-dtags), so a module
# LD_LIBRARY_PATH for the other architecture cannot override them.
# ==============================================================================
.DEFAULT_GOAL := all
LAMMPS_DIR ?= $(HOME)/lammps_latest
LAMMPSINC = $(LAMMPS_DIR)/fortran
LAMMPSLIB = $(LAMMPS_DIR)/build
NVHPC_ROOT ?= $(firstword $(wildcard /usr/local/modules/x64_v4/software/NVHPC/25.3-CUDA-12.8.0/Linux_x86_64/25.3 /opt/nvidia/hpc_sdk/Linux_x86_64/2025 /opt/nvidia/hpc_sdk/Linux_x86_64/25.9 /opt/nvidia/hpc_sdk/Linux_x86_64/current))
NVBIN ?= $(if $(EBROOTNVHPC),$(EBROOTNVHPC)/Linux_x86_64/$(EBVERSIONNVHPC)/compilers/bin,$(if $(NVHPC_ROOT),$(NVHPC_ROOT)/compilers/bin))
NVINCLUDE ?= $(if $(EBROOTNVHPC),$(EBROOTNVHPC)/Linux_x86_64/$(EBVERSIONNVHPC)/compilers/include,$(if $(NVHPC_ROOT),$(NVHPC_ROOT)/compilers/include))
NVLIBS ?= $(if $(EBROOTNVHPC),$(EBROOTNVHPC)/Linux_x86_64/$(EBVERSIONNVHPC)/compilers/lib,$(if $(NVHPC_ROOT),$(NVHPC_ROOT)/compilers/lib))
ifneq ($(NVBIN),)
export PATH := $(NVBIN):$(PATH)
endif
ifeq ($(origin FC),default)
FC = nvfortran
endif
FC ?= nvfortran
BINDIR = ../bin
BUILDDIR = build/$(ARCH)
OBJNAMES = Definitions.o potential_functions.o Read_input_data_nml.o Tools.o \
           Subsweep_Energy_CUDA.o InitializeCheckerboard.o MCsweep_Checkerboard.o \
           Shift_Cells.o CB_system_energy_gpu.o energy.o Volume_move.o \
           UpdateCheckerboard.o netcdf_trajectory.o mod_avbmc.o mod_clusters.o \
           lammps_hmc_wrapper.o Main.o
OBJ = $(addprefix $(BUILDDIR)/,$(OBJNAMES))

# A matching stack is available on Ladon. On other systems retain the original
# environment/module overrides and allow both host builds to share the libraries.
CLUSTER_ROOT = /usr/local/modules/x64_$(ARCH)/software
CLUSTER_NF = $(CLUSTER_ROOT)/netCDF-Fortran/4.6.1-NVHPC-25.3-CUDA-12.8.0
CLUSTER_NC = $(CLUSTER_ROOT)/netCDF/4.9.2-NVHPC-25.3-CUDA-12.8.0
CLUSTER_H5 = $(CLUSTER_ROOT)/HDF5/1.14.3-NVHPC-25.3-CUDA-12.8.0
NETCDF_FORTRAN_ROOT ?= $(if $(wildcard $(CLUSTER_NF)/include),$(CLUSTER_NF),$(if $(EBROOTNETCDFMINFORTRAN),$(EBROOTNETCDFMINFORTRAN),$(firstword $(wildcard /usr/local/netcdf-nv))))
NETCDF_C_ROOT ?= $(if $(wildcard $(CLUSTER_NC)/lib64),$(CLUSTER_NC),$(if $(EBROOTNETCDF),$(EBROOTNETCDF),$(NETCDF_FORTRAN_ROOT)))
NETCDFINC ?= $(if $(NETCDF_FORTRAN_ROOT),$(NETCDF_FORTRAN_ROOT)/include)
NETCDFLIB ?= $(firstword $(wildcard $(NETCDF_FORTRAN_ROOT)/lib64 $(NETCDF_FORTRAN_ROOT)/lib))
NETCDF_C_LIB ?= $(firstword $(wildcard $(NETCDF_C_ROOT)/lib64 $(NETCDF_C_ROOT)/lib))
HDF5LIB ?= $(firstword $(wildcard $(CLUSTER_H5)/lib $(EBROOTHDF5)/lib $(EBROOTHDF5)/lib64))
# Replace inherited x64_v4 module paths for v3, preserving explicit overrides.
ifeq ($(ARCH),v3)
ifneq ($(origin NETCDFINC),command line)
ifneq ($(filter /usr/local/modules/x64_v4/software/%,$(NETCDFINC)),)
NETCDFINC := $(CLUSTER_NF)/include
endif
endif
ifneq ($(origin NETCDFLIB),command line)
ifneq ($(filter /usr/local/modules/x64_v4/software/%,$(NETCDFLIB)),)
NETCDFLIB := $(CLUSTER_NF)/lib64
endif
endif
ifneq ($(origin NETCDF_C_LIB),command line)
ifneq ($(filter /usr/local/modules/x64_v4/software/%,$(NETCDF_C_LIB)),)
NETCDF_C_LIB := $(CLUSTER_NC)/lib64
endif
endif
ifneq ($(origin HDF5LIB),command line)
ifneq ($(filter /usr/local/modules/x64_v4/software/%,$(HDF5LIB)),)
HDF5LIB := $(CLUSTER_H5)/lib
endif
endif
endif
FCINC = -I$(BUILDDIR) $(if $(NETCDFINC),-I$(NETCDFINC)) $(if $(NVINCLUDE),-I$(NVINCLUDE)) \
        $(if $(FFTWINC),-I$(FFTWINC)) -I$(LAMMPSINC)
FCOPTS = -O3 -tp=x86-64-$(ARCH) -gpu=cc75,cc86,cc89,cc90,cc120 -cudalib=curand -module $(BUILDDIR) $(FCINC)
comma := ,
RPATH_FLAGS = -Wl$(comma)--disable-new-dtags \
              -Wl$(comma)-rpath$(comma)$(LAMMPSLIB) \
              $(if $(NETCDFLIB),-Wl$(comma)-rpath$(comma)$(NETCDFLIB)) \
              $(if $(NETCDF_C_LIB),-Wl$(comma)-rpath$(comma)$(NETCDF_C_LIB)) \
              $(if $(HDF5LIB),-Wl$(comma)-rpath$(comma)$(HDF5LIB))
LKLIBS = $(if $(NETCDFLIB),-L$(NETCDFLIB)) $(if $(NETCDF_C_LIB),-L$(NETCDF_C_LIB)) \
         $(if $(HDF5LIB),-L$(HDF5LIB)) $(if $(NVLIBS),-L$(NVLIBS)) \
         $(if $(FFTWLIB),-L$(FFTWLIB)) -L$(LAMMPSLIB) -llammps -lnetcdff -lnetcdf
LKOPTS = -cuda -gpu=cc75,cc86,cc89,cc90,cc120 -cudalib=curand -tp=x86-64-$(ARCH) $(RPATH_FLAGS)

.PHONY: all v3 v4 build-bin clean test_swap
all: v3 v4 $(BINDIR)/MCCB-gpu $(BINDIR)/mc_gpu.exe
v3:
	$(MAKE) ARCH=v3 build-bin
v4:
	$(MAKE) ARCH=v4 build-bin
build-bin: $(BINDIR)/MCCB-gpu-$(ARCH)
$(BINDIR) $(BUILDDIR):
	mkdir -p $@
$(BINDIR)/MCCB-gpu-$(ARCH): $(OBJ) $(BUILDDIR)/lammps.o | $(BINDIR)
	$(FC) $(LKOPTS) $(EXE_FLAG) -o $@ $(OBJ) $(BUILDDIR)/lammps.o $(LKLIBS)
$(BINDIR)/MCCB-gpu: MCCB-gpu.launcher | $(BINDIR)
	cp $< $@
	chmod +x $@
$(BINDIR)/mc_gpu.exe: $(BINDIR)/MCCB-gpu
	ln -sf MCCB-gpu $@
$(BUILDDIR)/%.o: %.cuf | $(BUILDDIR)
	$(FC) -c $(FCOPTS) $< -o $@
$(BUILDDIR)/lammps.o: $(LAMMPSINC)/lammps.f90 | $(BUILDDIR)
	$(FC) -c $(FCOPTS) $< -o $@

$(BUILDDIR)/potential_functions.o: $(BUILDDIR)/Definitions.o
$(BUILDDIR)/Subsweep_Energy_CUDA.o: $(BUILDDIR)/Definitions.o $(BUILDDIR)/potential_functions.o
$(BUILDDIR)/Tools.o: $(BUILDDIR)/Definitions.o $(BUILDDIR)/potential_functions.o
$(BUILDDIR)/Read_input_data_nml.o: $(BUILDDIR)/Definitions.o $(BUILDDIR)/Tools.o
$(BUILDDIR)/InitializeCheckerboard.o $(BUILDDIR)/MCsweep_Checkerboard.o: $(BUILDDIR)/Definitions.o $(BUILDDIR)/Subsweep_Energy_CUDA.o $(BUILDDIR)/Tools.o
$(BUILDDIR)/Shift_Cells.o $(BUILDDIR)/CB_system_energy_gpu.o $(BUILDDIR)/UpdateCheckerboard.o: $(BUILDDIR)/Definitions.o $(BUILDDIR)/Subsweep_Energy_CUDA.o
$(BUILDDIR)/energy.o: $(BUILDDIR)/Definitions.o $(BUILDDIR)/potential_functions.o $(BUILDDIR)/Tools.o
$(BUILDDIR)/Volume_move.o: $(BUILDDIR)/Definitions.o $(BUILDDIR)/Subsweep_Energy_CUDA.o $(BUILDDIR)/Tools.o
$(BUILDDIR)/netcdf_trajectory.o $(BUILDDIR)/mod_avbmc.o: $(BUILDDIR)/Definitions.o
$(BUILDDIR)/mod_clusters.o: $(BUILDDIR)/Definitions.o $(BUILDDIR)/mod_avbmc.o $(BUILDDIR)/Subsweep_Energy_CUDA.o
$(BUILDDIR)/lammps_hmc_wrapper.o: $(BUILDDIR)/Definitions.o $(BUILDDIR)/lammps.o
$(BUILDDIR)/Main.o: $(filter-out $(BUILDDIR)/Main.o,$(OBJ)) $(BUILDDIR)/lammps.o
clean:
	rm -rf build
	rm -f *.o *.mod *.smod test_swap_energy.exe $(BINDIR)/MCCB-gpu-v3 $(BINDIR)/MCCB-gpu-v4 $(BINDIR)/MCCB-gpu $(BINDIR)/mc_gpu.exe
# Legacy single-architecture test target.
test_swap: v3
	$(MAKE) ARCH=v3 test-swap-v3
test-swap-v3:
	$(FC) $(FCOPTS) $(LKOPTS) $(BUILDDIR)/Definitions.o $(BUILDDIR)/potential_functions.o $(BUILDDIR)/energy.o test_swap_energy.cuf -o test_swap_energy.exe
	./test_swap_energy.exe
	rm -f test_swap_energy.exe
