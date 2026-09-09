#!/bin/bash
# Master launcher for NPT vs NVT Consistency Benchmarks
# Executes simulations concurrently across ladon28 and ladon29

BASE_DIR="/home/e.lomba/MC_Checkerboard/test_cases/npt_nvt_consistency"

echo "Submitting Hard Spheres consistency job to ladon28..."
cd $BASE_DIR/hs && sbatch run_hs.slurm

echo "Submitting Lennard-Jones consistency job to ladon28..."
cd $BASE_DIR/lj && sbatch run_lj.slurm

echo "Submitting Site-Site Patchy consistency job to ladon29..."
cd $BASE_DIR/ssp && sbatch run_ssp.slurm

echo ""
echo "All consistency jobs submitted! Monitor with: squeue -u \$USER"
echo "When completed, evaluate consistency with: python3 $BASE_DIR/analyze_consistency.py $BASE_DIR"
