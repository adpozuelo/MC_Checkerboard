
# The following environment variables must be defined :
#
# NVBIN path of nvfortran binaries
# NVINCLUDE path to nvfortran header files and modules
# NVLIBS path to nvfortran dynamic libraries
# NETCDFINC path to NETCDF and NETCDFF includes
# NETCDFLIB path to NETCDF and NETCDFF libraries
# FFTWINC path to FFTW3 include files (might be in system's path)
# FFTWLIB path to FFTW2 dynamic libs (might be in system's path)
#
PATH := $(NVBIN):$(PATH)
FC = nvfortran
FCINC = -I$(NETCDFINC) -I$(NVINCLUDE) $(if $(FFTWINC),-I$(FFTWINC)) 
FCOPTS = -O3 -cudalib=curand $(FCINC)
LKLIBS = -L$(NETCDFLIB) -L$(NVLIBS) $(if $(FFTWLIB),-L$(FFTWLIB))
LKOPTS = -cuda -cudalib=curand  $(LKLIBS) -lnetcdff 

EXE = -o ../bin/mc_gpu.exe
OBJ = Definitions.o Read_input_data_nml.o Tools.o Subsweep_Energy_CUDA.o \
		InitializeCheckerboard.o MCsweep_Checkerboard.o Shift_Cells.o \
		CB_system_energy_gpu.o energy.o Volume_move.o UpdateCheckerboard.o \
		netcdf_trajectory.o \
		Main.o

%.o : %.mod
.SUFFIXES : .cuf .f90

%.o: %.cuf
	$(FC) -c $(FCOPTS) $< -o $@

all: $(OBJ) 
	$(FC) $(LKOPTS) $(LKINC) $(EXE) $(OBJ) $(LKLIBS)

clean:
	rm -f *.o *.mod ../bin/mc_gpu.exe
