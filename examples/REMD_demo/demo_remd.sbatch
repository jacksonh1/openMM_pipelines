#!/bin/bash
# Short end-to-end T-REMD demo: prepare() + remd.run() + analyze_remd() on a GPU.
#   cd examples/REMD_demo && sbatch demo_remd.sbatch     # submit from this folder
#SBATCH -p mit_normal_gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --gres=gpu:l40s:1
#SBATCH -t 00:30:00
#SBATCH -o demo_remd_%j.out
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate openmm_env

cd "$SLURM_SUBMIT_DIR"   # sbatch copies this script elsewhere; rely on the submit dir

echo "node: $(hostname)"
nvidia-smi -L
python ./demo_remd.py
