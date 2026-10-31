#!/bin/bash
# Resume/extend demo: a short md.run() continued to a longer target with resume().
#   cd examples/MD_demo && sbatch resume_demo.sbatch     # submit from this folder
# Walltime-recovery pattern: resubmitting this SAME script after a kill continues the run
# (resume_demo.py's run_or_resume checks for resume_manifest.json).
#SBATCH -p mit_normal_gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --gres=gpu:l40s:1
#SBATCH -t 00:20:00
#SBATCH -o resume_demo_%j.out
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate openmm_env

cd "$SLURM_SUBMIT_DIR"   # sbatch copies this script elsewhere; rely on the submit dir

echo "node: $(hostname)"
nvidia-smi -L
python ./resume_demo.py
