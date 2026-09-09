#!/usr/bin/env python3
"""
Convert LAMMPS data file from reduced LJ units (sigma=1.0) to real units (sigma=3.0 A)
by scaling box boundaries and coordinates by factor 3.0.
"""
import sys

def convert(infile, outfile, factor=3.0):
    with open(infile, 'r') as f:
        lines = f.readlines()

    out_lines = []
    in_atoms = False
    for line in lines:
        stripped = line.strip()
        if "units = lj" in line:
            line = line.replace("units = lj", "units = real")
        elif "xlo xhi" in line or "ylo yhi" in line or "zlo zhi" in line:
            parts = stripped.split()
            lo = float(parts[0]) * factor
            hi = float(parts[1]) * factor
            label = " ".join(parts[2:])
            line = f"  {lo:.12f} {hi:.12f} {label}\n"
        elif stripped.startswith("Atoms"):
            in_atoms = True
            out_lines.append(line)
            continue
        elif in_atoms and stripped and not stripped.startswith("#"):
            parts = stripped.split()
            if len(parts) >= 5:
                atom_id = parts[0]
                atom_type = parts[1]
                x = float(parts[2]) * factor
                y = float(parts[3]) * factor
                z = float(parts[4]) * factor
                rest = " ".join(parts[5:])
                line = f"{atom_id} {atom_type} {x:.14f} {y:.14f} {z:.14f} {rest}\n"

        out_lines.append(line)

    with open(outfile, 'w') as f:
        f.writelines(out_lines)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: convert_lammps_lj_to_real.py <input.atoms> <output.atoms> [factor]")
        sys.exit(1)
    fac = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
    convert(sys.argv[1], sys.argv[2], fac)
