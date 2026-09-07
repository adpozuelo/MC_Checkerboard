#!/bin/bash
# ==============================================================================
# Master SLURM submission script for NAHS Swap Reversibility & Performance
# ==============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=============================================================================="
echo " Submitting NAHS Identity Swap Benchmark to SLURM Short Queues"
echo " Working directory: $SCRIPT_DIR"
echo " Date:              $(date)"
echo "=============================================================================="

# Submit job WITH identity swaps
cd "$SCRIPT_DIR/run_with_swap"
JOB_SWAP=$(sbatch submit.sh | awk '{print $NF}')
echo " [+] Submitted 'run_with_swap'    -> Job ID: $JOB_SWAP"

# Submit job WITHOUT identity swaps
cd "$SCRIPT_DIR/run_without_swap"
JOB_NOSWAP=$(sbatch submit.sh | awk '{print $NF}')
echo " [+] Submitted 'run_without_swap' -> Job ID: $JOB_NOSWAP"

cd "$SCRIPT_DIR"

echo "=============================================================================="
echo " Both jobs submitted successfully!"
echo " Monitor queue:      squeue -u $USER"
echo " When completed, run: python3 analyze_benchmark.py"
echo "=============================================================================="
