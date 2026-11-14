#!/usr/bin/env bash
#SBATCH --job-name=electroexo-wahab
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
if [[ ! -f examples/slurm/run_demo.sbatch ]]; then
    printf 'Submit this job from the ElectroExo checkout root.\n' >&2
    exit 2
fi

export ELECTROEXO_PYTHON="$PWD/examples/slurm/wahab_python.sh"
export ELECTROEXO_OUTPUT_ROOT="${ELECTROEXO_OUTPUT_ROOT:-/scratch/${USER:?}/electroexo}"
exec bash "$PWD/examples/slurm/run_demo.sbatch" "$@"
