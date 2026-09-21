#!/bin/bash
#SBATCH --job-name=sb_avbmc
#SBATCH --partition=gpu_v4_short
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:15:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

export OMP_NUM_THREADS=4
echo "=== Starting SB AVBMC 50:50 Simulation on $(hostname) at $(date) ==="
/home/e.lomba/MC_Checkerboard/bin/MCCB-gpu datos.nml > output.out 2>&1
exit_code=$?
echo "=== Finished at $(date) with exit code ${exit_code} ==="
exit ${exit_code}
