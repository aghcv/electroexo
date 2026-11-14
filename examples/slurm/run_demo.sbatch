#!/usr/bin/env bash
#SBATCH --job-name=electroexo-demo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
scenario="${1:-examples/scenario_demo.yml}"
output_root="${ELECTROEXO_OUTPUT_ROOT:-results/slurm}"
python_executable="${ELECTROEXO_PYTHON:-python}"
repo_dir="$PWD"

mkdir -p "$output_root"
output_root="$(cd "$output_root" && pwd)"
run_dir="$(mktemp -d "$output_root/run-${SLURM_JOB_ID:-local}.XXXXXX")"
cp -- "$scenario" "$run_dir/scenario.yml"
printf 'Run directory %s\n' "$run_dir"

export MPLBACKEND=Agg
export MPLCONFIGDIR="$run_dir/.matplotlib"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"

"$python_executable" "$repo_dir/examples/slurm/run_with_provenance.py" \
    "$run_dir/scenario.yml" "$run_dir" "$repo_dir"
