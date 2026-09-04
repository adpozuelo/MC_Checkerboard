# Tabulated Potential Binary Mixture Example

This example demonstrates the simulation of a simple binary fluid mixture interacting via LAMMPS potential tables on the GPU using checkerboard Monte Carlo.

## Source & Potential
- Converted from the LAMMPS binary mixture benchmark in `~/IPHS/liquid/table/`.
- Interaction: Mie 50-49 potentials with $N = 100,000$ points, style `RSQ`, $r_{\min} = 10^{-6}$, $r_{\max} = 1.2$.
- Pair potential files: `pot11.dat`, `pot12.dat`, `pot22.dat` (also symlinked as `mie11.dat`, `mie12.dat`, `mie22.dat`).
- Particles: $N = 10,976$ atoms ($7,318$ of type 1 and $3,658$ of type 2) in a cubic box of size $L = 23.9396629$.

## Execution
Run with:
```bash
mpirun -np 1 ../../bin/mc_gpu.exe datos.nml 0
```

## Validation
- Initial LAMMPS potential energy: $776.2968458$ ($0.07072675$ / particle).
- Initial MC_Checkerboard GPU energy: $776.3055686$ ($0.07072755$ / particle), matching LAMMPS to $< 0.001\%$.
