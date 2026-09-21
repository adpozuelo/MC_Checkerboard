#!/usr/bin/env python3
"""
replicate_restart.py
Replicates a LAMMPS ellipsoid restart data file into an (nx * ny * nz) supercell.
Preserves particle types, orientations (quaternions), and density.
"""

import sys
import os
import argparse
import numpy as np

def replicate_system(input_file, output_file, nx=3, ny=3, nz=3):
    print(f"Reading {input_file}...")
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Parse headers and sections
    natoms = 0
    nellipsoids = 0
    ntypes = 0
    xlo, xhi = 0.0, 0.0
    ylo, yhi = 0.0, 0.0
    zlo, zhi = 0.0, 0.0

    atoms_lines = []
    ellipsoids_lines = []
    reading_atoms = False
    reading_ellipsoids = False

    for line in lines:
        s = line.strip()
        if not s:
            continue
        if "atoms" in s and "types" not in s and "ellipsoids" not in s:
            natoms = int(s.split()[0])
        elif "ellipsoids" in s and "types" not in s:
            nellipsoids = int(s.split()[0])
        elif "atom types" in s:
            ntypes = int(s.split()[0])
        elif "xlo xhi" in s:
            parts = s.split()
            xlo, xhi = float(parts[0]), float(parts[1])
        elif "ylo yhi" in s:
            parts = s.split()
            ylo, yhi = float(parts[0]), float(parts[1])
        elif "zlo zhi" in s:
            parts = s.split()
            zlo, zhi = float(parts[0]), float(parts[1])
        elif s.startswith("Atoms"):
            reading_atoms = True
            reading_ellipsoids = False
            continue
        elif s.startswith("Ellipsoids"):
            reading_atoms = False
            reading_ellipsoids = True
            continue
        elif reading_atoms:
            atoms_lines.append(s)
        elif reading_ellipsoids:
            ellipsoids_lines.append(s)

    Lx = xhi - xlo
    Ly = yhi - ylo
    Lz = zhi - zlo

    print(f"Original system:")
    print(f"  Atoms: {natoms}, Ellipsoids: {nellipsoids}, Types: {ntypes}")
    print(f"  Box: Lx={Lx:.6f}, Ly={Ly:.6f}, Lz={Lz:.6f}")
    print(f"  x in [{xlo:.6f}, {xhi:.6f}], y in [{ylo:.6f}, {yhi:.6f}], z in [{zlo:.6f}, {zhi:.6f}]")

    if len(atoms_lines) != natoms:
        raise ValueError(f"Expected {natoms} atoms, found {len(atoms_lines)}")
    if len(ellipsoids_lines) != nellipsoids:
        raise ValueError(f"Expected {nellipsoids} ellipsoids, found {len(ellipsoids_lines)}")

    # Parse atom coordinates
    # Format: id ityp eflag den x y z
    atoms_data = []
    for line in atoms_lines:
        parts = line.split()
        aid = int(parts[0])
        ityp = int(parts[1])
        eflag = int(parts[2])
        den = float(parts[3])
        x = float(parts[4])
        y = float(parts[5])
        z = float(parts[6])
        atoms_data.append((aid, ityp, eflag, den, x, y, z))

    # Sort by ID
    atoms_data.sort(key=lambda a: a[0])

    # Parse ellipsoids
    # Format: id sx sy sz qw qi qj qk
    ellips_data = []
    for line in ellipsoids_lines:
        parts = line.split()
        aid = int(parts[0])
        sx, sy, sz = float(parts[1]), float(parts[2]), float(parts[3])
        qw, qi, qj, qk = float(parts[4]), float(parts[5]), float(parts[6]), float(parts[7])
        ellips_data.append((aid, sx, sy, sz, qw, qi, qj, qk))

    ellips_data.sort(key=lambda e: e[0])

    # New system specifications
    total_replicas = nx * ny * nz
    new_natoms = natoms * total_replicas
    new_Lx = nx * Lx
    new_Ly = ny * Ly
    new_Lz = nz * Lz
    new_xlo = -new_Lx / 2.0
    new_xhi =  new_Lx / 2.0
    new_ylo = -new_Ly / 2.0
    new_yhi =  new_Ly / 2.0
    new_zlo = -new_Lz / 2.0
    new_zhi =  new_Lz / 2.0

    print(f"\nReplicated ({nx}x{ny}x{nz}) system:")
    print(f"  Total Replicas: {total_replicas}")
    print(f"  New Atoms: {new_natoms}")
    print(f"  New Box: Lx={new_Lx:.6f}, Ly={new_Ly:.6f}, Lz={new_Lz:.6f}")
    print(f"  New Bounds: x in [{new_xlo:.6f}, {new_xhi:.6f}], y in [{new_ylo:.6f}, {new_yhi:.6f}], z in [{new_zlo:.6f}, {new_zhi:.6f}]")

    # Shift offsets centered at 0
    x_shifts = [(ix - (nx - 1) / 2.0) * Lx for ix in range(nx)]
    y_shifts = [(iy - (ny - 1) / 2.0) * Ly for iy in range(ny)]
    z_shifts = [(iz - (nz - 1) / 2.0) * Lz for iz in range(nz)]

    print(f"  x shifts: {x_shifts}")
    print(f"  y shifts: {y_shifts}")
    print(f"  z shifts: {z_shifts}")

    new_atoms = []
    new_ellips = []

    atom_counter = 0
    for ix in range(nx):
        sx_shift = x_shifts[ix]
        for iy in range(ny):
            sy_shift = y_shifts[iy]
            for iz in range(nz):
                sz_shift = z_shifts[iz]
                for i in range(natoms):
                    atom_counter += 1
                    _, ityp, eflag, den, x, y, z = atoms_data[i]
                    _, sx, sy, sz, qw, qi, qj, qk = ellips_data[i]

                    new_x = x + sx_shift
                    new_y = y + sy_shift
                    new_z = z + sz_shift

                    new_atoms.append((atom_counter, ityp, eflag, den, new_x, new_y, new_z))
                    new_ellips.append((atom_counter, sx, sy, sz, qw, qi, qj, qk))

    assert len(new_atoms) == new_natoms
    assert len(new_ellips) == new_natoms

    # Write output LAMMPS data file
    print(f"\nWriting {output_file}...")
    with open(output_file, "w") as f:
        f.write("LAMMPS data file generated by replicate_restart.py (3x3x3 supercell)\n\n")
        f.write(f" {new_natoms:8d} atoms\n")
        f.write(f" {new_natoms:8d} ellipsoids\n")
        f.write(f" {ntypes:8d} atom types\n\n")
        f.write(f" {new_xlo:15.6f} {new_xhi:15.6f} xlo xhi\n")
        f.write(f" {new_ylo:15.6f} {new_yhi:15.6f} ylo yhi\n")
        f.write(f" {new_zlo:15.6f} {new_zhi:15.6f} zlo zhi\n\n")
        f.write("Atoms # ellipsoid\n\n")
        for a in new_atoms:
            f.write(f"{a[0]:8d} {a[1]:8d} {a[2]:3d} {a[3]:9.4f} {a[4]:16.8f} {a[5]:16.8f} {a[6]:16.8f}\n")
        f.write("\nEllipsoids\n\n")
        for e in new_ellips:
            f.write(f"{e[0]:8d} {e[1]:9.4f} {e[2]:9.4f} {e[3]:9.4f} {e[4]:16.8f} {e[5]:16.8f} {e[6]:16.8f} {e[7]:16.8f}\n")

    print(f"Successfully wrote {output_file} with {new_natoms} atoms.")

def main():
    parser = argparse.ArgumentParser(description="Replicate LAMMPS ellipsoid configuration")
    parser.add_argument("--input", "-i", default="../sb_avbmc_50_50/data.restart", help="Input restart file")
    parser.add_argument("--output", "-o", default="data.atoms", help="Output data file")
    parser.add_argument("--nx", type=int, default=3, help="Replication in x")
    parser.add_argument("--ny", type=int, default=3, help="Replication in y")
    parser.add_argument("--nz", type=int, default=3, help="Replication in z")

    args = parser.parse_args()
    replicate_system(args.input, args.output, args.nx, args.ny, args.nz)

if __name__ == "__main__":
    main()
