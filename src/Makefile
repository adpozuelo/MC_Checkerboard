PATH := $(NVBIN):$(PATH)
FC = nvfortran
FCINC = -I$(NETCDFINC) -I$(NVINCLUDE) $(if $(FFTWINC),-I$(FFTWINC)) 
FCOPTS = -O3 -cudalib=curand $(FCINC)
LKLIBS = -L$(NETCDFLIB) -L$(NVLIBS) $(if $(FFTWLIB),-L$(FFTWLIB))
LKOPTS = -cuda -cudalib=curand  $(LKLIBS)-lnetcdff -lnetcdf

EXE = -o ../bin/mc_gpu.exe
OBJ = Definitions.o Read_input_data.o Subsweep_Energy_CUDA.o \
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
