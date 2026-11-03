#!/bin/bash
#SBATCH --job-name=single_cell_tree_ilp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=5:00:00          # Gurobi limit 4 h plus model building and the heuristic start
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=abronner@princeton.edu
#SBATCH --output=/scratch/gpfs/MONA/abronner/slurm_logs/%j_log.txt
#SBATCH --error=/scratch/gpfs/MONA/abronner/slurm_logs/%j_error.txt

module load anaconda3/2024.6
module load gurobi/8.0.1
conda activate torch_cobra

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"             # must be the project root

EXP="${1:-single_cell_tree_ilp}"   # config name in configs_for_experiments/inference
echo "EXP=$EXP"
echo "EXP=$EXP" >&2
python -u single_cell_inference/snapshot_tree_ilp.py --config configs_for_experiments/inference/${EXP}.txt
