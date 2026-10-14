#!/bin/bash
#SBATCH --job-name=twotier_avbmc
#SBATCH --partition=gpu_v4_short
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:15:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

export OMP_NUM_THREADS=4
echo "=== Starting Two-Tier AVBMC test on $(hostname) at $(date) ==="
/home/e.lomba/MC_Checkerboard/bin/MCCB-gpu datos.nml > output.out 2>&1
echo "=== Finished at $(date) with exit code $? ==="
