FC = nvfortran
FCOPTS = -O3 -cudalib=curand
LKOPTS = -cuda -cudalib=curand

EXE = -o ../bin/mc_gpu.exe
OBJ = Definitions.o Read_input_data.o Subsweep_Energy_CUDA.o \
		InitializeCheckerboard.o MCsweep_Checkerboard.o Shift_Cells.o \
		CB_system_energy_gpu.o energy.o Volume_move.o UpdateCheckerboard.o \
		Main.o

%.o : %.mod
.SUFFIXES : .cuf .f90

%.o: %.cuf
	$(FC) -c $(FCOPTS) $< -o $@

all: $(OBJ) 
	$(FC) $(LKOPTS) $(LKINC) $(EXE) $(OBJ) $(LKLIBS)

clean:
	rm -f *.o *.mod ../bin/mc_gpu.exe
