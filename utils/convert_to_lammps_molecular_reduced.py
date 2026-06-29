#!/usr/bin/env python3
import sys
import os
import math

def main():
    if len(sys.argv) < 3:
        print("Usage: python convert_to_lammps_molecular_reduced.py <input.atoms> <output.atoms_lammps>")
        print("Outputs coordinates directly in reduced/LJ units (molecular style, compatible with mc_gpu.exe)")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    
    sigma = 1.0

    if not os.path.exists(input_path):
        print(f"Error: file '{input_path}' not found.")
        sys.exit(1)

    # Tetrahedral unit vectors
    s = 1.0 / math.sqrt(3.0)
    tetra_vecs = [
        [s, s, s],
        [s, -s, -s],
        [-s, s, -s],
        [-s, -s, s]
    ]

    print(f"Reading configuration from '{input_path}'...")
    with open(input_path, 'r') as f:
        lines = f.readlines()

    atoms_dict = {}
    ellipsoids_dict = {}
    xlo, xhi = 0.0, 0.0
    ylo, yhi = 0.0, 0.0
    zlo, zhi = 0.0, 0.0
    num_atoms = 0

    in_atoms = False
    in_ellipsoids = False

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        
        if 'atoms' in stripped and 'atom types' not in stripped:
            num_atoms = int(stripped.split()[0])
            continue
        if 'xlo xhi' in stripped:
            parts = stripped.split()
            xlo, xhi = float(parts[0]), float(parts[1])
            continue
        if 'ylo yhi' in stripped:
            parts = stripped.split()
            ylo, yhi = float(parts[0]), float(parts[1])
            continue
        if 'zlo zhi' in stripped:
            parts = stripped.split()
            zlo, zhi = float(parts[0]), float(parts[1])
            continue
            
        if stripped.startswith('Atoms'):
            in_atoms = True
            in_ellipsoids = False
            continue
        if stripped.startswith('Ellipsoids'):
            in_atoms = False
            in_ellipsoids = True
            continue
            
        if in_atoms:
            parts = stripped.split()
            atom_id = int(parts[0])
            itype = int(parts[1])
            x = float(parts[4])
            y = float(parts[5])
            z = float(parts[6])
            atoms_dict[atom_id] = (itype, x, y, z)
            
        if in_ellipsoids:
            parts = stripped.split()
            atom_id = int(parts[0])
            qw = float(parts[4])
            qi = float(parts[5])
            qj = float(parts[6])
            qk = float(parts[7])
            ellipsoids_dict[atom_id] = (qw, qi, qj, qk)

    Lx = xhi - xlo
    Ly = yhi - ylo
    Lz = zhi - zlo

    print(f"Writing molecular-style configuration (reduced/LJ units) to '{output_path}'...")
    with open(output_path, 'w') as f:
        f.write("# LAMMPS molecular data file for site-site patchy model (5-site representation, reduced/LJ units, centered)\n\n")
        f.write(f"{num_atoms * 5} atoms\n")
        f.write("4 atom types\n\n")
        
        f.write(f"{xlo:.8f} {xhi:.8f} xlo xhi\n")
        f.write(f"{ylo:.8f} {yhi:.8f} ylo yhi\n")
        f.write(f"{zlo:.8f} {zhi:.8f} zlo zhi\n\n")
        
        f.write("Atoms # molecular\n\n")
        
        new_atom_id = 1
        for atom_id in sorted(atoms_dict.keys()):
            itype, x_red, y_red, z_red = atoms_dict[atom_id]
            qw, qi, qj, qk = ellipsoids_dict[atom_id]
            
            center_type = 1 if itype == 1 else 3
            patch_type = 2 if itype == 1 else 4
            
            cx = x_red - Lx * math.floor(x_red / Lx + 0.5)
            cy = y_red - Ly * math.floor(y_red / Ly + 0.5)
            cz = z_red - Lz * math.floor(z_red / Lz + 0.5)
            
            # center site
            f.write(f"{new_atom_id} {atom_id} {center_type} {cx:.8f} {cy:.8f} {cz:.8f}\n")
            new_atom_id += 1
            
            qw2, qi2, qj2, qk2 = qw*qw, qi*qi, qj*qj, qk*qk
            qwi, qwj, qwk = qw*qi, qw*qj, qw*qk
            qij, qik, qjk = qi*qj, qi*qk, qj*qk
            
            r00 = qw2 + qi2 - qj2 - qk2
            r01 = 2.0 * (qij - qwk)
            r02 = 2.0 * (qik + qwj)
            
            r10 = 2.0 * (qij + qwk)
            r11 = qw2 - qi2 + qj2 - qk2
            r12 = 2.0 * (qjk - qwi)
            
            r20 = 2.0 * (qik - qwj)
            r21 = 2.0 * (qjk + qwi)
            r22 = qw2 - qi2 - qj2 + qk2
            
            for u in tetra_vecs:
                rx = r00 * u[0] + r01 * u[1] + r02 * u[2]
                ry = r10 * u[0] + r11 * u[1] + r12 * u[2]
                rz = r20 * u[0] + r21 * u[1] + r22 * u[2]
                
                px = cx + 0.5 * rx
                py = cy + 0.5 * ry
                pz = cz + 0.5 * rz
                
                # patch site
                f.write(f"{new_atom_id} {atom_id} {patch_type} {px:.8f} {py:.8f} {pz:.8f}\n")
                new_atom_id += 1

    print(f"Success! Wrote {num_atoms * 5} atoms to {output_path}.")

if __name__ == "__main__":
    main()
