#!/bin/bash
#SBATCH --job-name=consist_lj
#SBATCH --partition=gpu_v4_short
#SBATCH --nodelist=ladon28
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=slurm_lj_%j.out
#SBATCH --error=slurm_lj_%j.err

source /home/e.lomba/bin/setup_trj
export OMP_NUM_THREADS=4
BIN="/home/e.lomba/MC_Checkerboard/bin/MCCB-gpu"
BASE_DIR="/home/e.lomba/MC_Checkerboard/test_cases/npt_nvt_consistency/lj"

echo "=== [1/3] Running Lennard-Jones NpT Simulation on $(hostname) at $(date) ==="
cd $BASE_DIR/npt
$BIN datos.nml 0 > output.out 2>&1

echo "=== [2/3] Setting up NVT simulation from equilibrated NpT restart ==="
cd $BASE_DIR/nvt
cp $BASE_DIR/npt/data.restart data.atoms
cp $BASE_DIR/datos_nvt.nml datos.nml

echo "=== [3/3] Running Lennard-Jones NVT Simulation on $(hostname) at $(date) ==="
$BIN datos.nml 0 > output.out 2>&1

echo "=== Finished Lennard-Jones Consistency Run at $(date) ==="
