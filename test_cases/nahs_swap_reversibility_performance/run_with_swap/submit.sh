#!/bin/bash
#SBATCH --job-name=NAHS_with_swap
#SBATCH --partition=gpu_v4_short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH -o output_%j.out
#SBATCH -e error_%j.err

echo "=========================================================="
echo " Starting NAHS simulation WITH identity swaps"
echo " Job ID: $SLURM_JOB_ID on $HOSTNAME"
echo " Date:   $(date)"
echo "=========================================================="

source /usr/local/etc/modules.sh
module purge
source ~/bin/setup_trj

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
BIN="/home/e.lomba/MC_Checkerboard/bin/MCCB-gpu"

$BIN datos.nml

echo "=========================================================="
echo " Simulation completed at: $(date)"
echo "=========================================================="
