#!/bin/bash
#SBATCH --job-name=inference_errsummary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=abronner@princeton.edu
#SBATCH --output=/scratch/gpfs/MONA/abronner/slurm_logs/%j_errsummary_log.txt
#SBATCH --error=/scratch/gpfs/MONA/abronner/slurm_logs/%j_errsummary_error.txt

# End-of-array finalize job: aggregates the per-task error records each array task left under
# <manifest>.errparts into an errors-only summary written into each targeted experiment directory
# (inferred_models/<experiment>). It also compares the manifest against the on-disk output to list
# incomplete-but-unrecorded problems (tasks killed/timed out before finishing), querying sacct on the array
# job id for why those tasks ended. Submitted by submit_inference_array.sh with an afterany dependency on the
# array job, so it runs once every array task has finished (errored or not).
# Usage: sbatch summarize_inference_errors.sbatch <manifest> <chunk_size> <array_job_id>

module load anaconda3/2024.6
module load gurobi/8.0.1
conda activate torch_cobra

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

MANIFEST="${1:?Usage: sbatch summarize_inference_errors.sbatch <manifest> <chunk_size> <array_job_id>}"
CHUNK_SIZE="${2:?Usage: sbatch summarize_inference_errors.sbatch <manifest> <chunk_size> <array_job_id>}"
ARRAY_JID="${3:?Usage: sbatch summarize_inference_errors.sbatch <manifest> <chunk_size> <array_job_id>}"
BASE="$(basename "${MANIFEST%.*}")"
# bare filename: written into each experiment dir (inferred_models/<experiment>) by --summarize_errors
SUMMARY="errors_summary_${BASE#manifest_}.txt"

python -m inference.run_inference_on_data \
    --run_manifest "$MANIFEST" \
    --summarize_errors "$SUMMARY" \
    --chunk_size "$CHUNK_SIZE" \
    --array_job_id "$ARRAY_JID"
