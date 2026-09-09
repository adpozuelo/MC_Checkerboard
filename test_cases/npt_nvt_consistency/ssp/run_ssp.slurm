#!/bin/bash
#SBATCH --job-name=consist_ssp
#SBATCH --partition=gpu_v4_short
#SBATCH --nodelist=ladon29
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=slurm_ssp_%j.out
#SBATCH --error=slurm_ssp_%j.err

source /home/e.lomba/bin/setup_trj
export OMP_NUM_THREADS=4
BIN="/home/e.lomba/MC_Checkerboard/bin/MCCB-gpu"
BASE_DIR="/home/e.lomba/MC_Checkerboard/test_cases/npt_nvt_consistency/ssp"

echo "=== [1/3] Running Site-Site Patchy NpT Simulation on $(hostname) at $(date) ==="
cd $BASE_DIR/npt
$BIN datos.nml 0 > output.out 2>&1

echo "=== [2/3] Setting up NVT simulation from equilibrated NpT restart ==="
cd $BASE_DIR/nvt
cp $BASE_DIR/npt/data.restart data.atoms
cp $BASE_DIR/datos_nvt.nml datos.nml

echo "=== [3/3] Running Site-Site Patchy NVT Simulation on $(hostname) at $(date) ==="
$BIN datos.nml 0 > output.out 2>&1

echo "=== Finished Site-Site Patchy Consistency Run at $(date) ==="
