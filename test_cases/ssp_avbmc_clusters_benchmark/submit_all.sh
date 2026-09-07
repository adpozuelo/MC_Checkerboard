#!/bin/bash
# Master launcher for comparative AVBMC benchmarks
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Submitting WITH AVBMC benchmark..."
cd "$DIR/run_with_avbmc" && sbatch run.slurm

echo "Submitting WITHOUT AVBMC benchmark..."
cd "$DIR/run_without_avbmc" && sbatch run.slurm

echo "Both jobs submitted. Current queue status:"
squeue -u $USER
