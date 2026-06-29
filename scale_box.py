#!/usr/bin/env python3
import sys
import os

def main():
    # Print the warning immediately to highlight constraints
    print("=" * 80, file=sys.stderr)
    print(" WARNING: This script scales box bounds and atom coordinates relative to the origin.", file=sys.stderr)
    print("          This is ONLY valid for atomic or center-of-mass/rigid-body configurations.", file=sys.stderr)
    print("          If your system contains explicit bonds, angles, or molecular topology,", file=sys.stderr)
    print("          this scaling WILL modify bond lengths and structural constraints!", file=sys.stderr)
    print("=" * 80 + "\n", file=sys.stderr)

    if len(sys.argv) < 3:
        print("Usage: python scale_box.py <filename.atoms> <scaling_factor>")
        sys.exit(1)

    input_file = sys.argv[1]
    try:
        factor = float(sys.argv[2])
    except ValueError:
        print("Error: scaling_factor must be a float value.")
        sys.exit(1)

    if not os.path.exists(input_file):
        print(f"Error: file '{input_file}' not found.")
        sys.exit(1)

    base, ext = os.path.splitext(input_file)
    output_file = f"{base}_scaled{ext}"

    print(f"Scaling box and coordinates in '{input_file}' by a factor of {factor}...")
    print(f"Writing output to '{output_file}'...")

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        in_atoms = False
        coord_indices = None

        for line in infile:
            stripped = line.strip()

            # Handle blank lines or headers
            if not stripped:
                outfile.write(line)
                continue

            # Detect Atoms section
            if stripped.startswith("Atoms"):
                in_atoms = True
                outfile.write(line)
                continue
            elif in_atoms and stripped and stripped[0].isalpha():
                # End of Atoms section
                in_atoms = False
                outfile.write(line)
                continue

            # Process box bounds
            if not in_atoms and any(keyword in stripped for keyword in ["xlo xhi", "ylo yhi", "zlo zhi"]):
                parts = stripped.split()
                if len(parts) >= 4:
                    try:
                        lo = float(parts[0]) * factor
                        hi = float(parts[1]) * factor
                        outfile.write(f" {lo:.8f} {hi:.8f} {parts[2]} {parts[3]}\n")
                        continue
                    except ValueError:
                        pass
                outfile.write(line)
                continue

            # Process atoms coordinates
            if in_atoms:
                parts = stripped.split()
                if len(parts) >= 5:
                    if coord_indices is None:
                        # Determine coordinate columns based on field count
                        n_fields = len(parts)
                        if n_fields == 5:
                            coord_indices = (2, 3, 4)
                        elif n_fields == 6:
                            coord_indices = (3, 4, 5)
                        else:
                            # 7 or more fields
                            coord_indices = (4, 5, 6)
                    
                    try:
                        idx_x, idx_y, idx_z = coord_indices
                        parts[idx_x] = f"{float(parts[idx_x]) * factor:.8f}"
                        parts[idx_y] = f"{float(parts[idx_y]) * factor:.8f}"
                        parts[idx_z] = f"{float(parts[idx_z]) * factor:.8f}"
                        # Write formatted line
                        outfile.write(" " + " ".join(parts) + "\n")
                        continue
                    except (ValueError, IndexError):
                        pass
                outfile.write(line)
                continue

            outfile.write(line)

    print("Success! Scale transformation complete.")

if __name__ == "__main__":
    main()
