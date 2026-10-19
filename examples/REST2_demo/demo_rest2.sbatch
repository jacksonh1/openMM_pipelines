#!/bin/bash
# Short end-to-end REST2 demo: prepare() + rest2.run() + analyze_rest2() on a GPU.
#   cd examples/REST2_demo && sbatch demo_rest2.sbatch   # submit from this folder
#SBATCH -p mit_normal_gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --gres=gpu:l40s:1
#SBATCH -t 00:30:00
#SBATCH -o demo_rest2_%j.out
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate openmm_env

cd "$SLURM_SUBMIT_DIR"   # sbatch copies this script elsewhere; rely on the submit dir

echo "node: $(hostname)"
nvidia-smi -L
python ./demo_rest2.py
