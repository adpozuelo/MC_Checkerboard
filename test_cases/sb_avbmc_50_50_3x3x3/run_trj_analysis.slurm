#!/bin/bash
#SBATCH --job-name=trj_3x3x3
#SBATCH --partition=gpu_v4_short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --output=trj_analysis.out

set -euo pipefail
source /usr/local/modules/x64_v4/software/Lmod/8.7.59/lmod/lmod/init/bash
module use /usr/local/modules/x64_v4/modules/all
module load netCDF/4.9.2-NVHPC-25.3-CUDA-12.8.0
module load netCDF-Fortran/4.6.1-NVHPC-25.3-CUDA-12.8.0

cd /home/e.lomba/MC_Checkerboard/test_cases/sb_avbmc_50_50_3x3x3
/sif/users/e.lomba/trj_analysis/bin/trj_analysis trj.nml 0
