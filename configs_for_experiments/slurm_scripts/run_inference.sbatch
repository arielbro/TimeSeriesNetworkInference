#!/bin/bash
#SBATCH --job-name=model_inference
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=256G # the noise robustness experiment is the bottleneck
#SBATCH --time=5:00:00
#SBATCH --output=slurm_logs/slurm-%j.out
#SBATCH --error=slurm_logs/slurm-%j.err
#SBATCH --mail-type=FAIL,END        # send email when job begins
#SBATCH --mail-user=abronner@princeton.edu

#SBATCH --output=/scratch/gpfs/MONA/abronner/slurm_logs/%j_log.txt
#SBATCH --error=/scratch/gpfs/MONA/abronner/slurm_logs/%j_error.txt

module load anaconda3/2024.6
module load gurobi/8.0.1
conda activate torch_cobra

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"             # must be the project root

# export OMP_NUM_THREADS=1
# export OPENBLAS_NUM_THREADS=1
# export MKL_NUM_THREADS=1
# Make Gurobi honor the allocation
EXP="${1:?Usage: sbatch run_experiment.sbatch <experiment_name>}"  # change per experiment
echo "EXP=$EXP"
echo "EXP=$EXP" >&2

# Step 2: run inference on it
python -m inference.run_inference_on_data \
    --config configs_for_experiments/inference/${EXP}.txt
