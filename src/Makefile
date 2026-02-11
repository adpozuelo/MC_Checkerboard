FC = nvfortran
FCOPTS = -O3 -cudalib=curand -I/usr/local/modules/x64_v4/software/netCDF-Fortran/4.6.1-NVHPC-25.3-CUDA-12.8.0/include
LKOPTS = -cuda -cudalib=curand -L/usr/local/modules/x64_v4/software/netCDF-Fortran/4.6.1-NVHPC-25.3-CUDA-12.8.0/lib -lnetcdff -L/usr/local/modules/x64_v4/software/netCDF/4.9.2-NVHPC-25.3-CUDA-12.8.0/lib64 -lnetcdf

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
