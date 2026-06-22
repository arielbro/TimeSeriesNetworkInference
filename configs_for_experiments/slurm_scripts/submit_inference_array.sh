#!/bin/bash
# Enumerate every single-graph inference problem for an inference config (inference-parameter combinations x
# subexperiments x graphs), then submit a SLURM array with one task per chunk of problems.
#
# Run from the project root:
#   configs_for_experiments/slurm_scripts/submit_inference_array.sh configs_for_experiments/inference/<exp>.txt

# ---- per-run resources (edit here; passed to sbatch, overriding the run_inference_array.sbatch defaults) ----
CPUS=4                 # CPUs per task, fed to Gurobi as its Threads parameter
TIME=04:00:00          # wall time per task; must cover (chunk_size x model_inference_timeout_secs) + overhead
MEM=16G                # memory per task
MAX_CONCURRENT=200     # cap on simultaneously running array tasks (the %K throttle)
MAX_ARRAY_SIZE=2501    # cluster MaxArraySize (check: scontrol show config | grep MaxArraySize)
# ------------------------------------------------------------------------------------------------------------

CONFIG="${1:?Usage: submit_inference_array.sh <inference_config.txt>}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

module load anaconda3/2024.6
module load gurobi/8.0.1
conda activate torch_cobra
set -euo pipefail

# manifest name keyed to the config so concurrent runs of different configs don't collide
MANIFEST="manifest_$(basename "${CONFIG%.*}").jsonl"
echo "Enumerating problems for ${CONFIG} ..."
python -m inference.run_inference_on_data --config "$CONFIG" --emit_manifest "$MANIFEST"
N=$(wc -l < "$MANIFEST")
echo "Enumerated ${N} problems -> ${MANIFEST}"
if [ "$N" -eq 0 ]; then
    echo "No problems to run; not submitting."
    exit 0
fi

# choose the smallest chunk size that keeps the number of array tasks within MaxArraySize
CHUNK_SIZE=$(( (N + MAX_ARRAY_SIZE - 1) / MAX_ARRAY_SIZE ))
N_TASKS=$(( (N + CHUNK_SIZE - 1) / CHUNK_SIZE ))
echo "chunk_size=${CHUNK_SIZE}  n_tasks=${N_TASKS}  max_concurrent=${MAX_CONCURRENT}"

sbatch --array=0-$((N_TASKS - 1))%${MAX_CONCURRENT} \
       --cpus-per-task="${CPUS}" --time="${TIME}" --mem="${MEM}" \
       "${SCRIPT_DIR}/run_inference_array.sbatch" "$MANIFEST" "$CHUNK_SIZE"
