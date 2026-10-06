#!/bin/bash
# Short end-to-end demo of prepare() + md.run() on a GPU.
#   cd examples/MD_demo && sbatch demo_short_run.sbatch     # submit from this folder
#SBATCH -p mit_normal_gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --gres=gpu:l40s:1
#SBATCH -t 00:20:00
#SBATCH -o demo_short_run_%j.out
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate openmm_env

cd "$SLURM_SUBMIT_DIR"   # sbatch copies this script elsewhere; rely on the submit dir

echo "node: $(hostname)"
nvidia-smi -L
python ./demo_short_run.py
