#!/usr/bin/env bash
#SBATCH --job-name=avbmc_graph_fix_zero_eps
#SBATCH --partition=gpu_v4_short
#SBATCH --gres=gpu:nvidia_rtx_pro_4500_blackwell:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
sha256sum --check /sif/users/e.lomba/MC_Checkerboard/bin/avbmc_graph_fix_probe/binary.sha256
python3 case.py prepare
export OMP_NUM_THREADS=2
export MCCB_LAUNCHER_DEBUG=1
export MCCB_AVBMC_GRAPH_CHECK=1
export MCCB_AVBMC_GRAPH_CHECK_FROM_STEP=0
/sif/users/e.lomba/MC_Checkerboard/bin/avbmc_graph_fix_probe/MCCB-gpu-v4 datos.nml > simulation.out 2>&1
python3 case.py analyze > analysis.out
