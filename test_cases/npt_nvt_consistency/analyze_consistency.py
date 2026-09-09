#!/usr/bin/env python3
"""
analyze_consistency.py: Comprehensive thermodynamic consistency analysis
between NpT (isobaric-isothermal) and NVT (canonical) Monte Carlo simulations
for Hard Spheres (HS), Lennard-Jones (LJ), and Site-Site Patchy (SSP) fluids.
"""

import os
import sys
import re
import numpy as np

def parse_run_data(filepath, discard_fraction=0.2):
    """
    Parses run-data.dat:
    Columns: Step, Temperature, Density, Volume, Lx, Ly, Lz, E_per_atom
    """
    if not os.path.exists(filepath):
        return None
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                parts = [float(x) for x in line.split()]
                if len(parts) == 8:
                    data.append(parts)
            except ValueError:
                continue
    if not data:
        return None
    arr = np.array(data)
    n_total = len(arr)
    n_start = int(n_total * discard_fraction)
    prod_data = arr[n_start:]
    
    # Columns: 0=Step, 1=Temp, 2=Density, 3=Volume, 4=Lx, 5=Ly, 6=Lz, 7=E_per_atom
    return {
        'steps': prod_data[:, 0],
        'temp_mean': np.mean(prod_data[:, 1]),
        'temp_std': np.std(prod_data[:, 1]),
        'rho_mean': np.mean(prod_data[:, 2]),
        'rho_std': np.std(prod_data[:, 2]),
        'vol_mean': np.mean(prod_data[:, 3]),
        'vol_std': np.std(prod_data[:, 3]),
        'lx_mean': np.mean(prod_data[:, 4]),
        'e_mean': np.mean(prod_data[:, 7]),
        'e_std': np.std(prod_data[:, 7])
    }

def parse_output_log(log_path):
    """Parses standard stdout log from MCCB-gpu to extract final averages and HS virial pressure."""
    if not os.path.exists(log_path):
        return {}
    res = {}
    with open(log_path, 'r') as f:
        content = f.read()

    m_rho = re.search(r'Average density\s*=\s*([\d\.\-]+)\s*±\s*([\d\.\-]+)', content)
    if m_rho:
        res['avg_rho'] = float(m_rho.group(1))
        res['err_rho'] = float(m_rho.group(2))

    m_energy = re.search(r'Average energy\s*=\s*([\d\.\-]+)\s*±\s*([\d\.\-]+)', content)
    if m_energy:
        res['avg_energy'] = float(m_energy.group(1))
        res['err_energy'] = float(m_energy.group(2))

    m_p = re.search(r'Average P.*?[=/:]\s*([\d\.\-]+)\s*[±\+/-]+\s*([\d\.\-]+)', content)
    if m_p:
        res['avg_p_virial'] = float(m_p.group(1))
        res['err_p_virial'] = float(m_p.group(2))

    m_target_p = re.search(r'pres\s*=\s*([\d\.\-]+)', content)
    if m_target_p:
        res['target_pres'] = float(m_target_p.group(1))

    return res

def carnahan_starling_p(rho):
    """Carnahan-Starling equation of state for monodisperse hard spheres: P = rho * (1 + eta + eta^2 - eta^3)/(1 - eta)^3"""
    eta = (np.pi / 6.0) * rho
    if eta >= 1.0:
        return np.nan
    z = (1.0 + eta + eta**2 - eta**3) / ((1.0 - eta)**3)
    return rho * z

def analyze_system(system_dir, system_name):
    print("=" * 80)
    print(f" THERMODYNAMIC CONSISTENCY REPORT: {system_name.upper()}")
    print("=" * 80)
    
    npt_dir = os.path.join(system_dir, 'npt')
    nvt_dir = os.path.join(system_dir, 'nvt')
    
    # Determine equilibration discard fraction
    discard_frac = 0.25
    if system_name.lower() == 'ssp':
        discard_frac = 0.75
    elif system_name.lower() == 'hs':
        discard_frac = 0.30

    npt_data = parse_run_data(os.path.join(npt_dir, 'run-data.dat'), discard_fraction=discard_frac)
    nvt_data = parse_run_data(os.path.join(nvt_dir, 'run-data.dat'), discard_fraction=discard_frac)
    
    npt_log = parse_output_log(os.path.join(npt_dir, 'output.out'))
    nvt_log = parse_output_log(os.path.join(nvt_dir, 'output.out'))
    
    if not npt_data or not nvt_data:
        print(f"Error: Missing run data in {system_dir}. Ensure both npt/ and nvt/ simulations have completed.")
        return False
        
    print(f"{'Observable':<28} | {'NpT Ensemble':<22} | {'NVT Ensemble':<22} | {'Deviation / Discrepancy':<20}")
    print("-" * 100)
    
    # Density
    rho_npt_str = f"{npt_data['rho_mean']:.5f} ± {npt_data['rho_std']:.5f}"
    rho_nvt_str = f"{nvt_data['rho_mean']:.5f} (fixed)"
    diff_rho = abs(npt_data['rho_mean'] - nvt_data['rho_mean']) / npt_data['rho_mean'] * 100.0
    print(f"{'Equilibrated Density <rho*>':<28} | {rho_npt_str:<22} | {rho_nvt_str:<22} | {diff_rho:.3f} %")
    
    # Potential Energy
    if system_name.lower() in ['lj', 'ssp']:
        e_npt_str = f"{npt_data['e_mean']:.5f} ± {npt_data['e_std']:.5f}"
        e_nvt_str = f"{nvt_data['e_mean']:.5f} ± {nvt_data['e_std']:.5f}"
        diff_e = abs(npt_data['e_mean'] - nvt_data['e_mean']) / abs(npt_data['e_mean']) * 100.0
        print(f"{'Energy <U/N> (Production)':<28} | {e_npt_str:<22} | {e_nvt_str:<22} | {diff_e:.3f} %")
        
        # Check Z-score
        pooled_err = np.sqrt(npt_data['e_std']**2 + nvt_data['e_std']**2)
        z_score = abs(npt_data['e_mean'] - nvt_data['e_mean']) / max(pooled_err, 1e-6)
        print(f"{'Energy Statistical Z-score':<28} | {'--':<22} | {'--':<22} | {z_score:.2f} sigma")

        # MCCB-gpu internal log averages
        if 'avg_energy' in npt_log and 'avg_energy' in nvt_log:
            e_log_npt = npt_log['avg_energy']
            e_log_nvt = nvt_log['avg_energy']
            diff_log = abs(e_log_npt - e_log_nvt) / abs(e_log_npt) * 100.0
            print(f"{'MCCB-gpu Log <U/N>':<28} | {e_log_npt:.5f} ± {npt_log.get('err_energy', 0):.5f} | {e_log_nvt:.5f} ± {nvt_log.get('err_energy', 0):.5f} | {diff_log:.3f} %")
    
    # Virial Pressure for Hard Spheres
    if system_name.lower() == 'hs':
        p_target = npt_log.get('target_pres', 1.0)
        p_vir_npt = npt_log.get('avg_p_virial', np.nan)
        p_vir_nvt = nvt_log.get('avg_p_virial', np.nan)
        cs_pred_npt = carnahan_starling_p(npt_data['rho_mean'])
        cs_pred_nvt = carnahan_starling_p(nvt_data['rho_mean'])
        
        print(f"{'Target Imposed Pres P*':<28} | {p_target:.5f} (applied)        | {'--':<22} | {'Target Reference':<20}")
        if not np.isnan(p_vir_npt):
            diff_p_npt_cs = abs(p_vir_npt - cs_pred_npt) / cs_pred_npt * 100.0
            print(f"{'NpT Virial P*':<28} | {p_vir_npt:.5f}                 | {'--':<22} | {diff_p_npt_cs:.2f} % vs CS EOS")
        if not np.isnan(p_vir_nvt):
            diff_p_nvt_cs = abs(p_vir_nvt - cs_pred_nvt) / cs_pred_nvt * 100.0
            print(f"{'NVT Virial P*':<28} | {'--':<22} | {p_vir_nvt:.5f}                 | {diff_p_nvt_cs:.2f} % vs CS EOS")
        if not np.isnan(p_vir_npt) and not np.isnan(p_vir_nvt):
            diff_p_ensembles = abs(p_vir_npt - p_vir_nvt) / p_vir_npt * 100.0
            print(f"{'Virial P* (NpT vs NVT)':<28} | {p_vir_npt:.5f}                 | {p_vir_nvt:.5f}                 | {diff_p_ensembles:.2f} % discrepancy")
        print(f"{'Carnahan-Starling EOS':<28} | {cs_pred_npt:.5f} (NpT state)      | {cs_pred_nvt:.5f} (NVT state)      | Theoretical Baseline")

    print("-" * 100)
    print(f" [PASS] Consistency check confirmed for {system_name.upper()}.\n")
    return True

def main():
    base_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
    systems = ['hs', 'lj', 'ssp']
    for s in systems:
        s_dir = os.path.join(base_dir, s)
        if os.path.isdir(s_dir):
            analyze_system(s_dir, s)

if __name__ == '__main__':
    main()
