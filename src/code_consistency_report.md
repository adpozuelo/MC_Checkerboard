# Code Consistency & Bug Fix Report

This report documents the consistency issues and bugs identified in the codebase and records the fixes that have been applied.

---

## 1. GPU Memory Crash during Volume Moves (NpT Ensemble) — **FIXED**
### Location
* File: [Subsweep_Energy_CUDA.cuf](file:///Users/elomba/MC_Checkerboard/src/Subsweep_Energy_CUDA.cuf) (subroutines `clean`, `energyinitialize`, and `subsweepgpu`)
* File: [Volume_move.cuf](file:///Users/elomba/MC_Checkerboard/src/Volume_move.cuf)

### Description
In NpT simulations, rebuilding the cell grid via `energyinitialize` calls `clean()`, which previously:
1. Deallocated constant patchy potential arrays on the GPU (`Npatches_dev`, `Ntor_angle_dev`, `patchdev`, `ref_tor_dev`, `tor_angle_dev`, `Vpot_Matrix_dev`, `epsilon_LJ_dev`), causing illegal GPU memory access errors on subsequent energy evaluations.
2. Failed to reset `init_subsweep` and `init_shift` saved module flags, leading to skipped allocations of array parameters like `Offset_dev` and subsequent crashes when assigning to them.

### Applied Fix
* Modified `clean()` to retain constant potential parameters on the GPU.
* Reset `init_subsweep = 0` and `init_shift = 0` in `clean()` to trigger safe re-allocation of all cell-dependent arrays on the next Monte Carlo sweep / cell boundary shift.

---

## 2. NetCDF Write Variable Mismatch (Uninitialized ID) — **FIXED**
### Location
* File: [netcdf_trajectory.cuf](file:///Users/elomba/MC_Checkerboard/src/netcdf_trajectory.cuf)

### Description
The definition of the per-atom stress variable (`c_stress`) was commented out inside `netcdf_init_trajectory`, leaving `c_stress_varid` uninitialized, whereas `netcdf_write_frame` actively attempted to write stress data using `c_stress_varid`.

### Applied Fix
* Uncommented the definition of `c_stress` in `netcdf_init_trajectory` to register the variable correctly and match the schema expected during trajectory write operations.

---

## 3. Syntax Error in `energy.cuf` — **FIXED**
### Location
* File: [energy.cuf](file:///Users/elomba/MC_Checkerboard/src/energy.cuf#L316)

### Description
Line 316 contained a stray trailing parenthesis (`Else if (pot_int == 1) then)`) which caused compile-time syntax errors.

### Applied Fix
* Removed the stray parenthesis.
