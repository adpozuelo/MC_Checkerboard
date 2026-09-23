#!/usr/bin/env python3
"""
Setup script for Sanchez-Burgos 50:50 Scaffold + Surfactant AVBMC simulation.
Extracts tabulated potentials from potentials.table and converts LAMMPS configuration.
"""
import os
import sys
import numpy as np

def rot_to_quat(R):
    tr = R[0, 0] + R[1, 1] + R[2, 2]
    if tr > 0:
        S = np.sqrt(tr + 1.0) * 2.0
        qw = 0.25 * S
        qx = (R[2, 1] - R[1, 2]) / S
        qy = (R[0, 2] - R[2, 0]) / S
        qz = (R[1, 0] - R[0, 1]) / S
    elif (R[0, 0] > R[1, 1]) and (R[0, 0] > R[2, 2]):
        S = np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2.0
        qw = (R[2, 1] - R[1, 2]) / S
        qx = 0.25 * S
        qy = (R[0, 1] + R[1, 0]) / S
        qz = (R[0, 2] + R[2, 0]) / S
    elif R[1, 1] > R[2, 2]:
        S = np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2.0
        qw = (R[0, 2] - R[2, 0]) / S
        qx = (R[0, 1] + R[1, 0]) / S
        qy = 0.25 * S
        qz = (R[1, 2] + R[2, 1]) / S
    else:
        S = np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2.0
        qw = (R[1, 0] - R[0, 1]) / S
        qx = (R[0, 2] + R[2, 0]) / S
        qy = (R[1, 2] + R[2, 1]) / S
        qz = 0.25 * S
    q = np.array([qw, qx, qy, qz])
    q = q / np.linalg.norm(q)
    if q[0] < 0:
        q = -q
    return q

def extract_tables(table_file, out_dir):
    print(f"Reading {table_file}...")
    with open(table_file, "r") as f:
        lines = f.readlines()
        
    sections = {}
    current_sec = None
    current_lines = []
    
    for line in lines:
        s = line.strip()
        if s in ["PHS", "CSW", "ZERO"]:
            if current_sec is not None:
                sections[current_sec] = current_lines
            current_sec = s
            current_lines = [line]
        elif current_sec is not None:
            current_lines.append(line)
            
    if current_sec is not None:
        sections[current_sec] = current_lines
        
    print(f"Found sections: {list(sections.keys())}")
    
    # Mapping for MC_Checkerboard pot_int == 3 (alternating type indices: 1=core0, 2=patch0, 3=core1, 4=patch1)
    # pot11.dat: Core 0 - Core 0 (PHS)
    # pot13.dat: Core 0 - Core 1 (PHS)
    # pot33.dat: Core 1 - Core 1 (PHS)
    # pot24.dat: Patch 0 - Patch 1 (CSW, unlike patches)
    # pot22.dat: Patch 0 - Patch 0 (ZERO, like patches)
    # pot44.dat: Patch 1 - Patch 1 (ZERO, like patches)
    # pot12, pot14, pot23, pot34: Core-Patch (ZERO)
    
    table_map = {
        "pot11.dat": "PHS",
        "pot13.dat": "PHS",
        "pot33.dat": "PHS",
        "pot24.dat": "CSW",
        "pot22.dat": "CSW",
        "pot44.dat": "ZERO",
        "pot12.dat": "ZERO",
        "pot14.dat": "ZERO",
        "pot23.dat": "ZERO",
        "pot34.dat": "ZERO",
    }
    
    for fname, sec in table_map.items():
        fpath = os.path.join(out_dir, fname)
        with open(fpath, "w") as f:
            f.write(f"# Table file for pair {fname} ({sec})\n")
            f.writelines(sections[sec])
        print(f"  Wrote {fpath} ({sec}, {len(sections[sec])} lines)")

def convert_configuration(lammps_data_path, out_dir):
    print(f"Reading LAMMPS configuration: {lammps_data_path}")
    with open(lammps_data_path, "r") as f:
        lines = f.readlines()
        
    xlo, xhi = 0.0, 0.0
    ylo, yhi = 0.0, 0.0
    zlo, zhi = 0.0, 0.0
    atom_lines = []
    reading_atoms = False
    
    for line in lines:
        parts = line.strip().split()
        if not parts:
            continue
        if "xlo" in line and "xhi" in line:
            xlo, xhi = float(parts[0]), float(parts[1])
        elif "ylo" in line and "yhi" in line:
            ylo, yhi = float(parts[0]), float(parts[1])
        elif "zlo" in line and "zhi" in line:
            zlo, zhi = float(parts[0]), float(parts[1])
        elif parts[0] == "Atoms":
            reading_atoms = True
            continue
        elif parts[0] == "Velocities":
            reading_atoms = False
            continue
        elif reading_atoms:
            atom_lines.append(line)
            
    Lx = xhi - xlo
    Ly = yhi - ylo
    Lz = zhi - zlo
    L = np.array([Lx, Ly, Lz])
    print(f"Box dimensions: Lx={Lx:.8f}, Ly={Ly:.8f}, Lz={Lz:.8f}")
    
    atoms_by_mol = {}
    for line in atom_lines:
        parts = line.strip().split()
        if len(parts) >= 6:
            aid = int(parts[0])
            mid = int(parts[1])
            atype = int(parts[2])
            pos = np.array([float(parts[3]), float(parts[4]), float(parts[5])])
            if mid not in atoms_by_mol:
                atoms_by_mol[mid] = []
            atoms_by_mol[mid].append((aid, atype, pos))
            
    num_mols = len(atoms_by_mol)
    print(f"Total molecules: {num_mols}")
    
    inv_sq3 = 1.0 / np.sqrt(3.0)
    scaffold_ref = np.array([
        [+inv_sq3, +inv_sq3, +inv_sq3],
        [+inv_sq3, -inv_sq3, -inv_sq3],
        [-inv_sq3, +inv_sq3, -inv_sq3],
        [-inv_sq3, -inv_sq3, +inv_sq3]
    ])
    
    surfactant_ref = np.array([
        [1.0, 0.0, 0.0],
        [-0.5, np.sqrt(3.0)/2.0, 0.0],
        [-0.5, -np.sqrt(3.0)/2.0, 0.0]
    ])
    
    molecules = []
    max_fit_err = 0.0
    
    for mid in sorted(atoms_by_mol.keys()):
        matoms = sorted(atoms_by_mol[mid], key=lambda a: a[0])
        core_list = [a for a in matoms if a[1] in [1, 2]]
        patch_list = [a for a in matoms if a[1] in [3, 4]]
        
        assert len(core_list) == 1, f"Molecule {mid} does not have exactly 1 core"
        core = core_list[0]
        patches = patch_list
        core_pos = core[2]
        core_type = core[1]
        
        if core_type == 1:
            mc_type = 0 # Scaffold
            ref_dirs = scaffold_ref
            n_patches = 4
        elif core_type == 2:
            mc_type = 1 # Surfactant
            ref_dirs = surfactant_ref
            n_patches = 3
        else:
            raise ValueError(f"Unknown core type: {core_type}")
            
        assert len(patches) == n_patches, f"Molecule {mid}: expected {n_patches} patches, got {len(patches)}"
        
        actual_dirs = []
        for p in patches:
            dr = p[2] - core_pos
            dr = dr - L * np.round(dr / L)
            r_patch = np.linalg.norm(dr)
            actual_dirs.append(dr / r_patch)
        actual_dirs = np.array(actual_dirs)
        
        H = ref_dirs.T @ actual_dirs
        U, S, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        D = np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
        
        fit_err = np.max(np.abs((ref_dirs @ R.T) - actual_dirs))
        if fit_err > max_fit_err:
            max_fit_err = fit_err
            
        q = rot_to_quat(R)
        
        core_wrapped = core_pos - np.array([xlo, ylo, zlo])
        core_wrapped = core_wrapped - L * np.floor(core_wrapped / L)
        
        molecules.append({
            'mid': mid,
            'mc_type': mc_type,
            'pos': core_wrapped,
            'quat': q
        })
        
    print(f"Maximum patch reconstruction error: {max_fit_err:.3e}")
    
    data_atoms_path = os.path.join(out_dir, "data.atoms")
    print(f"Writing {data_atoms_path}...")
    with open(data_atoms_path, "w") as f:
        f.write("LAMMPS data file generated for MC_Checkerboard\n\n")
        f.write(f"    {num_mols} atoms\n")
        f.write(f"    {num_mols} ellipsoids\n")
        f.write(f"       2 atom types\n\n")
        f.write(f"       0.00000000       {Lx:.8f} xlo xhi\n")
        f.write(f"       0.00000000       {Ly:.8f} ylo yhi\n")
        f.write(f"       0.00000000       {Lz:.8f} zlo zhi\n\n")
        
        f.write("Atoms # ellipsoid\n")
        for i, m in enumerate(molecules):
            aid = i + 1
            atype = m['mc_type'] + 1
            f.write(f"{aid:8d} {atype:6d} 1 1.0000 {m['pos'][0]:16.8f} {m['pos'][1]:16.8f} {m['pos'][2]:16.8f}\n")
            
        f.write("\nEllipsoids\n")
        for i, m in enumerate(molecules):
            aid = i + 1
            q = m['quat']
            f.write(f"{aid:8d} 1.0000 1.0000 1.0000 {q[0]:16.8f} {q[1]:16.8f} {q[2]:16.8f} {q[3]:16.8f}\n")
            
    print(f"Successfully wrote {data_atoms_path} with {num_mols} molecules.")

def write_datos_nml(out_dir):
    datos_path = os.path.join(out_dir, "datos.nml")
    print(f"Writing {datos_path}...")
    content = """&Control_Params
  istep_ini      = 0,
  istep_fin      = 50,
  Neq            = 0,
  Nmove          = 20,
  Nsave          = 2,
  Ndump          = 10,
  Nrestart       = 50,
  imovie         = .true.,
  traj_format    = 'netcdf',
  Npart_types    = 2,
  data_file      = 'data.atoms',
  Ecpu_check     = .true.,
  ncluster       = 5,
  rcl            = 1.20,
  minPts         = 3,
  cluster_types  = 1, 2,
  asym_threshold = 0.5,
  border_criterion = 'asymmetry',
  avbmc          = .true.,
  avbmc_k_trials = 10,
  avbmc_max_trials = 50
/
&MC_Params
  hmax         = 0.015,
  omax         = 0.1,
  displ_update = .true.,
  temp0        = 1.0,
  temp1        = 1.0,
  npt          = .false.,
  seed         = 492845
/
&Potential_Params
  model         = 'SB',
  sigma_jon_aux = -0.3,
  sigma_tor_jon = 0.0,
  rangepp       = 2.0,
  Bool_tor      = 0,
  patch_radial_factor = 0.5,
  eps_R         = 0.6666666666666667,
  eps_CSW       = 11.111111111111111,
  rw_CSW        = 0.12,
  alpha_CSW     = 0.005,
  rc_CSW        = 0.20,
  table_mc      = .true.,
  table_file_mc = 'pot'
/
--- TOPOLOGY AND PATCH GEOMETRY DATA ---
          4           3
  0.57735027  0.57735027  0.57735027  0.70710678 -0.70710678  0.00000000     0
  0.57735027 -0.57735027 -0.57735027  0.70710678  0.70710678  0.00000000     0
 -0.57735027  0.57735027 -0.57735027  0.70710678  0.70710678  0.00000000     0
 -0.57735027 -0.57735027  0.57735027  0.70710678 -0.70710678  0.00000000     0
  1.00000000  0.00000000  0.00000000  0.00000000  1.00000000  0.00000000     0
 -0.50000000  0.86602540  0.00000000  0.00000000  1.00000000  0.00000000     0
 -0.50000000 -0.86602540  0.00000000  0.00000000  1.00000000  0.00000000     0
-- vpot matrix --
 1.00 1.00 1.00 1.00 1.00 1.00 1.00
 1.00 1.00 1.00 1.00 1.00 1.00 1.00
 1.00 1.00 1.00 1.00 1.00 1.00 1.00
 1.00 1.00 1.00 1.00 1.00 1.00 1.00
 1.00 1.00 1.00 1.00 0.00 0.00 0.00
 1.00 1.00 1.00 1.00 0.00 0.00 0.00
 1.00 1.00 1.00 1.00 0.00 0.00 0.00
-- sigmas_ang --
 0.300 0.300 0.300 0.300
 0.300 0.300 0.300
-- sigma_LJ --
 1.000 1.000
 1.000 1.000
"""
    with open(datos_path, "w") as f:
        f.write(content)
    print(f"Successfully wrote {datos_path}.")

if __name__ == "__main__":
    out_dir = "/home/e.lomba/MC_Checkerboard/test_cases/sb_avbmc_50_50"
    table_src = "/home/e.lomba/sanchez_burgos_20260915/sanchez_burgos_lammps/potentials.table"
    config_src = "/home/e.lomba/sanchez_burgos_20260915/sanchez_burgos_lammps/cases/cube_N2000_one/system.data"
    
    if len(sys.argv) > 1:
        config_src = sys.argv[1]
        
    extract_tables(table_src, out_dir)
    convert_configuration(config_src, out_dir)
    write_datos_nml(out_dir)
    print("Setup completed successfully.")
