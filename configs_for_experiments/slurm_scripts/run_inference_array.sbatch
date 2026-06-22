#!/bin/bash
#SBATCH --job-name=inference_array
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=abronner@princeton.edu
#SBATCH --output=/scratch/gpfs/MONA/abronner/slurm_logs/%A_%a_log.txt
#SBATCH --error=/scratch/gpfs/MONA/abronner/slurm_logs/%A_%a_error.txt

# One SLURM array task runs a contiguous chunk of single-graph inference problems from a manifest.
# Resources (--cpus-per-task / --time / --mem / --array) are set by submit_inference_array.sh and override the
# #SBATCH defaults above. CPUs per task are fed to Gurobi as its Threads parameter (the inference code reads
# SLURM_CPUS_PER_TASK). Normally invoked via submit_inference_array.sh, not directly.
# Usage: sbatch run_inference_array.sbatch <manifest> <chunk_size>

module load anaconda3/2024.6
module load gurobi/8.0.1
conda activate torch_cobra

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# let the numeric libraries use the allocation; Gurobi's own thread count is set in code from SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

MANIFEST="${1:?Usage: sbatch run_inference_array.sbatch <manifest> <chunk_size>}"
CHUNK_SIZE="${2:?Usage: sbatch run_inference_array.sbatch <manifest> <chunk_size>}"

python -m inference.run_inference_on_data \
    --run_manifest "$MANIFEST" \
    --task_id "$SLURM_ARRAY_TASK_ID" \
    --chunk_size "$CHUNK_SIZE"
