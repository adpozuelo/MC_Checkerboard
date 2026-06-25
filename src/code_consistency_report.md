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

---

## 4. Incorrect Potential call in `accumulate_atom_energy_sitesite` — **FIXED**
### Location
* File: [Subsweep_Energy_CUDA.cuf](file:///Users/elomba/MC_Checkerboard/src/Subsweep_Energy_CUDA.cuf) (subroutine `accumulate_atom_energy_sitesite`)

### Description
For the site-site patchy model (`pot_int == 2`), the per-atom energy calculation subroutine `accumulate_atom_energy_sitesite` mistakenly called `vlj` (the standard Lennard-Jones core potential) instead of `vlj_palaia` (the WCA core + cosine tail potential). This caused a discrepancy between the pair/cell energy calculations (which correctly used `vlj_palaia`) and the per-atom energy calculation.

### Applied Fix
* Changed `call vlj` to `call vlj_palaia` inside `accumulate_atom_energy_sitesite`.

---

## 5. Colloid-Colloid Potential Sign Mismatch (`epsp`) in `vlj_palaia` — **FIXED**
### Location
* File: [Subsweep_Energy_CUDA.cuf](file:///Users/elomba/MC_Checkerboard/src/Subsweep_Energy_CUDA.cuf) (device function `vlj_palaia`)
* File: [energy.cuf](file:///Users/elomba/MC_Checkerboard/src/energy.cuf) (subroutine `vlj_palaia`)

### Description
`VL0` is calculated in `Read_input_data_nml.cuf` as `(sigma/rc)^12 - (sigma/rc)^6` (which evaluates to a negative number). This value is passed as `epsp` (the tail depth) to `vlj_palaia`. However, both the CPU and GPU implementations of `vlj_palaia` calculated the tail interaction as `-epsp * cos(...)**2` and the WCA minimum shift as `-epsp`. Because `epsp` was negative, the tail became positive (repulsive) and the minimum shift became positive.

### Applied Fix
* Modified `vlj_palaia` on both CPU and GPU to use the absolute value `abs(epsp)` for the tail depth, ensuring the potential remains physically correct (attractive tail and core shift downwards).

---

## 6. Single-Patch Configuration Parsing Bug — **FIXED**
### Location
* File: [Read_input_data_nml.cuf](file:///Users/elomba/MC_Checkerboard/src/Read_input_data_nml.cuf) (subroutine `Read_Input_Data`)

### Description
In systems where a particle type has exactly one patch (`Nrb_sites(J) == 1`), the reader loaded coordinates directly into `patch(Indx+K)` but then calculated `norm` using uninitialized/stale variables `px, py, pz` and wrote to a stale `Indx2` index. This corrupted the patch vectors and risked memory corruption.

### Applied Fix
* Rewrote the single-patch parsing logic to read directly into local variables `px, py, pz`, calculate `norm` correctly, and store the normalized values in the correct `Indx` location.
