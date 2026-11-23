#!/bin/bash
# Requeue auto-resume demo: submit ONCE; SLURM reruns the SAME job on preemption / node failure and
# run_or_resume() continues from the last checkpoint.
#   cd examples/MD_demo && sbatch requeue_demo.sbatch     # submit from this folder
# This is the whole requeue idiom: --requeue + a body that always calls run_or_resume. No
# self-resubmit, no SIGTERM trap. A non-zero crash is NOT requeued, so a run that cannot progress
# fails once (loudly) instead of looping.
#SBATCH -p mit_preemptable        # <- set to YOUR cluster's preemptable partition (the point of --requeue)
#SBATCH --requeue                 # let SLURM requeue this job on preemption / node failure
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --gres=gpu:l40s:1
#SBATCH -t 00:20:00
#SBATCH -o requeue_demo_%j.out
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate openmm_env

cd "$SLURM_SUBMIT_DIR"   # sbatch copies this script elsewhere; rely on the submit dir

echo "node: $(hostname)  requeue attempt: ${SLURM_RESTART_COUNT:-0}"
nvidia-smi -L
python ./requeue_demo.py
