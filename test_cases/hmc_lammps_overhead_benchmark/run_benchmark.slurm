#!/bin/bash
#SBATCH --job-name=hmc_bench
#SBATCH --partition=gpu_v4_short
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=slurm_bench_%j.out
#SBATCH --error=slurm_bench_%j.err

source /home/e.lomba/bin/setup_trj
export OPAL_PREFIX=/yggdrasil/usr_local/modules/x64_v4/software/NVHPC/25.3-CUDA-12.8.0/Linux_x86_64/25.3/comm_libs/12.8/hpcx/hpcx-2.22.1/ompi
export PATH=$OPAL_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$HOME/lammps_latest/build:$OPAL_PREFIX/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

echo "=== Starting Comprehensive HMC Overhead Benchmark on $(hostname) at $(date) ==="
$OPAL_PREFIX/bin/mpirun -np 1 ./benchmark_hmc_suite.exe datos.nml 0 > benchmark_output.out 2>&1
echo "=== Finished at $(date) with exit code $? ==="
