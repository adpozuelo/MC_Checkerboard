#!/usr/bin/env bash
#SBATCH --job-name=ssp_av40_graph_fix_replay
#SBATCH --partition=gpu_v4_long
#SBATCH --gres=gpu:nvidia_rtx_pro_4500_blackwell:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=logs/graph_fix_replay_%j.out
#SBATCH --error=logs/graph_fix_replay_%j.err
set -euo pipefail

campaign=${SLURM_SUBMIT_DIR:?Submit from growth_extension_2026_11}
sha256sum --check /sif/users/e.lomba/MC_Checkerboard/bin/avbmc_graph_fix_probe/binary.sha256
cd "$campaign/graph_fix_replay_1to1_av40_100k_seed356"
if [[ -e output.out ]]; then
  echo "Existing output.out; refusing to overwrite" >&2
  exit 3
fi
export OMP_NUM_THREADS=4
export MCCB_LAUNCHER_DEBUG=1
export MCCB_AVBMC_GRAPH_CHECK=1
export MCCB_AVBMC_GRAPH_CHECK_FROM_STEP=299260
/sif/users/e.lomba/MC_Checkerboard/bin/avbmc_graph_fix_probe/MCCB-gpu-v4 datos.nml > output.out 2>&1
grep -q 'Simulation completely finished' output.out
