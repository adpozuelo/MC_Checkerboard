#!/usr/bin/env python3
import sys
import os

def main():
    if len(sys.argv) < 3:
        print("Usage: python scale_trajectory.py <filename.lammpstrj> <sigma_in_A>")
        sys.exit(1)

    input_file = sys.argv[1]
    try:
        sigma = float(sys.argv[2])
    except ValueError:
        print("Error: sigma must be a float value.")
        sys.exit(1)

    if not os.path.exists(input_file):
        print(f"Error: file '{input_file}' not found.")
        sys.exit(1)

    # Determine output filename
    base, ext = os.path.splitext(input_file)
    output_file = f"{base}_real{ext}"

    print(f"Scaling '{input_file}' to real units (sigma = {sigma} A)...")
    print(f"Writing scaled trajectory to '{output_file}'...")

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        in_atoms = False
        in_box = 0
        atom_cols = []
        x_idx = y_idx = z_idx = -1

        for line in infile:
            stripped = line.strip()
            
            # Detect item headers
            if stripped.startswith("ITEM: BOX BOUNDS"):
                in_box = 3
                outfile.write(line)
                continue
            elif stripped.startswith("ITEM: ATOMS"):
                in_atoms = True
                in_box = 0
                # Parse columns to find coordinates x, y, z
                parts = stripped.split()
                atom_cols = parts[2:]
                try:
                    x_idx = atom_cols.index('x')
                    y_idx = atom_cols.index('y')
                    z_idx = atom_cols.index('z')
                except ValueError:
                    # Fallback to column index guess if 'x', 'y', 'z' are named differently
                    x_idx = next((i for i, c in enumerate(atom_cols) if 'x' in c.lower()), -1)
                    y_idx = next((i for i, c in enumerate(atom_cols) if 'y' in c.lower()), -1)
                    z_idx = next((i for i, c in enumerate(atom_cols) if 'z' in c.lower()), -1)
                
                outfile.write(line)
                continue
            elif stripped.startswith("ITEM:"):
                in_atoms = False
                in_box = 0
                outfile.write(line)
                continue

            # Process box bounds
            if in_box > 0:
                parts = stripped.split()
                if len(parts) >= 2:
                    try:
                        val_lo = float(parts[0]) * sigma
                        val_hi = float(parts[1]) * sigma
                        if len(parts) > 2: # Triclinic bounds
                            val_extra = float(parts[2]) * sigma
                            outfile.write(f"      {val_lo:.7f}     {val_hi:.7f}     {val_extra:.7f}\n")
                        else:
                            outfile.write(f"      {val_lo:.7f}     {val_hi:.7f}\n")
                    except ValueError:
                        outfile.write(line)
                else:
                    outfile.write(line)
                in_box -= 1
                continue

            # Process atom coordinates
            if in_atoms:
                parts = stripped.split()
                if len(parts) > max(x_idx, y_idx, z_idx) and x_idx != -1:
                    try:
                        parts[x_idx] = f"{float(parts[x_idx]) * sigma:.7f}"
                        parts[y_idx] = f"{float(parts[y_idx]) * sigma:.7f}"
                        parts[z_idx] = f"{float(parts[z_idx]) * sigma:.7f}"
                        outfile.write("        " + "      ".join(parts) + "\n")
                    except ValueError:
                        outfile.write(line)
                else:
                    outfile.write(line)
                continue

            outfile.write(line)

    print("Success! Scale transformation complete.")

if __name__ == "__main__":
    main()
