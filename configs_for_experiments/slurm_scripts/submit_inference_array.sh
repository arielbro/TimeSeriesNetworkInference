#!/bin/bash
# Enumerate every single-graph inference problem for an inference config (inference-parameter combinations x
# subexperiments x graphs), then submit a SLURM array with one task per chunk of problems.
#
# Run from the project root:
#   configs_for_experiments/slurm_scripts/submit_inference_array.sh configs_for_experiments/inference/<exp>.txt

# ---- per-run resources (edit here; passed to sbatch, overriding the run_inference_array.sbatch defaults) ----
CPUS=4                 # CPUs per task, fed to Gurobi as its Threads parameter
MEM=16G                # memory per task
TIME_OVERHEAD=1.05     # wall time = chunk_size x model_inference_timeout_secs x this (5% slack); TIME is computed below
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

# wall time per task = chunk_size x per-graph Gurobi timeout x TIME_OVERHEAD. A task runs CHUNK_SIZE problems
# serially, each capped at model_inference_timeout_secs; that timeout is constant across the config, so read it
# from the manifest the problems were enumerated into (taking the max guards against any per-problem variation).
TIME=$(python - "$MANIFEST" "$CHUNK_SIZE" "$TIME_OVERHEAD" <<'PY'
import json, math, sys
manifest, chunk_size, overhead = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
timeouts = set()
with open(manifest) as f:
    for line in f:
        if line.strip():
            timeouts.add(json.loads(line)["kwargs"].get("model_inference_timeout_secs"))
if not timeouts or None in timeouts:
    sys.exit("Cannot derive wall time: model_inference_timeout_secs is unset (no Gurobi time limit) in the manifest")
secs = max(60, int(math.ceil(max(timeouts) * chunk_size * overhead)))  # 60s floor for sane minimum
days, rem = divmod(secs, 86400)
h, rem = divmod(rem, 3600)
m, s = divmod(rem, 60)
print("{:d}-{:02d}:{:02d}:{:02d}".format(days, h, m, s) if days else "{:d}:{:02d}:{:02d}".format(h, m, s))
PY
)
echo "time=${TIME} (chunk_size x model_inference_timeout_secs x ${TIME_OVERHEAD})"

# clear any error records left by a previous submission of this config so the summary reflects only this run
rm -rf "${MANIFEST}.errparts"

ARRAY_JID=$(sbatch --parsable --array=0-$((N_TASKS - 1))%${MAX_CONCURRENT} \
       --cpus-per-task="${CPUS}" --time="${TIME}" --mem="${MEM}" \
       "${SCRIPT_DIR}/run_inference_array.sbatch" "$MANIFEST" "$CHUNK_SIZE")
echo "submitted inference array as job ${ARRAY_JID}"

# finalize job: runs once every array task has finished (afterany = regardless of success/failure) and writes
# the errors-only summary aggregating what the tasks recorded.
SUMMARY_JID=$(sbatch --parsable --dependency=afterany:"${ARRAY_JID}" \
       "${SCRIPT_DIR}/summarize_inference_errors.sbatch" "$MANIFEST" "$CHUNK_SIZE" "$ARRAY_JID")
echo "submitted error-summary finalize as job ${SUMMARY_JID} (afterany:${ARRAY_JID})"
