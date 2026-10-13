import numpy as np

def generate_lj_data(filename, N=1000, rho=0.03):
    V = N / rho
    L = V ** (1.0 / 3.0)
    half_L = L / 2.0
    
    # 10 x 10 x 10 lattice
    nx = ny = nz = 10
    dx = L / nx
    
    positions = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                x = -half_L + (ix + 0.5) * dx
                y = -half_L + (iy + 0.5) * dx
                z = -half_L + (iz + 0.5) * dx
                positions.append((x, y, z))
                
    positions = positions[:N]
    
    with open(filename, 'w') as f:
        f.write("LAMMPS data file for LJ fluid validation\n\n")
        f.write(f" {N:10d} atoms\n")
        f.write("          1 atom types\n\n")
        f.write(f" {-half_L:16.8f} {half_L:16.8f} xlo xhi\n")
        f.write(f" {-half_L:16.8f} {half_L:16.8f} ylo yhi\n")
        f.write(f" {-half_L:16.8f} {half_L:16.8f} zlo zhi\n\n")
        f.write("Atoms # atomic\n\n")
        for i, (x, y, z) in enumerate(positions, 1):
            f.write(f"{i:10d}         1 {x:16.8f} {y:16.8f} {z:16.8f}\n")

if __name__ == '__main__':
    generate_lj_data('data.atoms', 1000, 0.03)
    print("Generated data.atoms for N=1000, rho=0.03")
