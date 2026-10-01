#!/bin/bash
# Run the test suite on a GPU (the `slow` tests run real OpenMM dynamics, far faster
# on CUDA than CPU). Submit from the repo root so $SLURM_SUBMIT_DIR is the package:
#   sbatch -w node3620 tests/gpu_tests.sbatch        # or -w node3619
# Pass pytest args through PYTEST_ARGS:
#   sbatch -w node3620 --export=ALL,PYTEST_ARGS="tests/test_stages.py" tests/gpu_tests.sbatch
#SBATCH -p pi_keating
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --gres=gpu:l40s:1
#SBATCH -t 00:20:00
#SBATCH -o gpu_tests_%j.out
set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate openmm_env

# sbatch copies this script to a spool path, so BASH_SOURCE is useless here — rely on
# $SLURM_SUBMIT_DIR (the dir sbatch was called from) for the repo location.
cd "$SLURM_SUBMIT_DIR"

export OPENMM_TEST_PLATFORM=CUDA
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "node: $(hostname)  platform: $OPENMM_TEST_PLATFORM"
nvidia-smi -L
python -m openmm.testInstallation

python -m pytest -q ${PYTEST_ARGS:-}
