#!/usr/bin/env python3
"""
==============================================================================
Master Figure Reproduction Script for MCCB-gpu Manuscript
==============================================================================
Executes all figure plotting scripts in sequence and validates output files.

Manuscript Figures:
  - Figure 2: consistency_all.pdf (Thermodynamic Consistency across NpT and NVT)
  - Figure 3: scaling_benchmark.pdf (CPU vs GPU Scaling across Ensembles)
  - Figure 4: hmc_overhead_benchmark_comparison.pdf (In-Memory HMC Profiling)
  - Figure 5: avbmc_benchmark_comparison.pdf (AVBMC Growth Kinetics & Acceptance)
  - Figure 6: avbmc_csd_distribution.pdf (Final Cluster Size Distribution)

Usage:
  python reproduce_all_figures.py [--output-dir ./figures]

Author: Enrique Lomba / Antonio Diaz-Pozuelo / Eva G. Noya
==============================================================================
"""

import os
import sys
import argparse
import subprocess
import time

def check_dependencies():
    missing = []
    try:
        import numpy
    except ImportError:
        missing.append('numpy')
    try:
        import matplotlib
    except ImportError:
        missing.append('matplotlib')

    if missing:
        print(f"Error: Missing required Python dependencies: {', '.join(missing)}")
        print("Please install them using: pip install -r requirements.txt")
        sys.exit(1)

def run_script(script_path, args=None):
    cmd = [sys.executable, script_path]
    if args:
        cmd.extend(args)
    print(f"\n>> Running {os.path.basename(script_path)}...")
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    if res.returncode != 0:
        print(f"FAILED (elapsed: {dt:.2f}s):\n{res.stderr}")
        return False
    print(res.stdout.strip())
    print(f"Completed in {dt:.2f}s")
    return True

def main():
    check_dependencies()

    root_dir = os.path.dirname(os.path.abspath(__file__))
    scripts_dir = os.path.join(root_dir, 'scripts')
    data_dir = os.path.join(root_dir, 'data')
    default_out_dir = os.path.join(root_dir, 'figures')

    parser = argparse.ArgumentParser(description='Reproduce all figures for the MCCB-gpu manuscript')
    parser.add_argument('--output-dir', default=default_out_dir, help='Directory to store generated figures')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 80)
    print("  MCCB-GPU MANUSCRIPT: INDEPENDENT FIGURE REPRODUCTION SUITE")
    print("=" * 80)
    print(f"Root Directory   : {root_dir}")
    print(f"Data Directory   : {data_dir}")
    print(f"Output Directory : {args.output_dir}")

    tasks = [
        ('plot_fig2_consistency.py', [
            '--data-dir', os.path.join(data_dir, 'consistency'),
            '--out-dir', args.output_dir
        ]),
        ('plot_fig3_scaling.py', [
            '--data-dir', os.path.join(data_dir, 'scaling'),
            '--out-dir', args.output_dir
        ]),
        ('plot_fig4_hmc_overhead.py', [
            '--data-dir', os.path.join(data_dir, 'hmc'),
            '--out-dir', args.output_dir
        ]),
        ('plot_fig5_fig6_avbmc.py', [
            '--data-dir', os.path.join(data_dir, 'avbmc'),
            '--out-dir', args.output_dir
        ]),
    ]

    all_ok = True
    for script_name, script_args in tasks:
        script_path = os.path.join(scripts_dir, script_name)
        ok = run_script(script_path, script_args)
        if not ok:
            all_ok = False

    expected_files = [
        ('Fig 2 (all)', 'consistency_all.pdf'),
        ('Fig 2 (HS)',  'consistency_hs.pdf'),
        ('Fig 2 (LJ)',  'consistency_lj.pdf'),
        ('Fig 2 (SSP)', 'consistency_ssp.pdf'),
        ('Fig 3',       'scaling_benchmark.pdf'),
        ('Fig 4',       'hmc_overhead_benchmark_comparison.pdf'),
        ('Fig 5',       'avbmc_benchmark_comparison.pdf'),
        ('Fig 6',       'avbmc_csd_distribution.pdf'),
        ('AVBMC (all)', 'avbmc_benchmark_all.pdf'),
    ]

    print("\n" + "=" * 80)
    print("  REPRODUCTION SUMMARY & ARTIFACT VERIFICATION")
    print("=" * 80)
    print(f"{'Manuscript Reference':<20} | {'Expected File':<36} | {'Status':<10} | {'Size'}")
    print("-" * 80)

    for ref, fname in expected_files:
        fpath = os.path.join(args.output_dir, fname)
        png_path = os.path.join(args.output_dir, fname.replace('.pdf', '.png'))
        if os.path.exists(fpath) and os.path.getsize(fpath) > 0:
            sz_pdf = f"{os.path.getsize(fpath) / 1024:.1f} KB"
            status = "OK (PDF+PNG)" if os.path.exists(png_path) else "OK (PDF)"
            print(f"{ref:<20} | {fname:<36} | {status:<10} | {sz_pdf}")
        else:
            print(f"{ref:<20} | {fname:<36} | MISSING    | -")
            all_ok = False

    print("=" * 80)
    if all_ok:
        print("[SUCCESS] All manuscript figures were reproduced successfully!")
    else:
        print("[ERROR] One or more figures could not be reproduced. See error log above.")
        sys.exit(1)

if __name__ == '__main__':
    main()
