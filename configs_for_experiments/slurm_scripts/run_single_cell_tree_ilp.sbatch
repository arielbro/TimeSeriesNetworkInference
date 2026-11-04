#!/bin/bash
#SBATCH --job-name=single_cell_tree_ilp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=5:00:00          # Gurobi limit 4 h plus model building and the heuristic start
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=abronner@princeton.edu
#SBATCH --output=/scratch/gpfs/MONA/abronner/slurm_logs/%A_%a_log.txt
#SBATCH --error=/scratch/gpfs/MONA/abronner/slurm_logs/%A_%a_error.txt

# Runs single_cell_inference/snapshot_tree_ilp.py on a config from configs_for_experiments/inference. As a
# SLURM array (submit_single_cell_tree_ilp.sh), task i runs only the i-th method of the config's method list;
# otherwise every method runs in turn.
# Usage: sbatch [--array=0-<n_methods - 1>] run_single_cell_tree_ilp.sbatch [<config name>]

module load anaconda3/2024.6
module load gurobi/8.0.1
conda activate torch_cobra

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"             # must be the project root

EXP="${1:-single_cell_tree_ilp}"
CONFIG="configs_for_experiments/inference/${EXP}.txt"
echo "EXP=$EXP"
echo "EXP=$EXP" >&2

METHOD_ARGS=()
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    METHOD=$(python -c "import ast, sys; v = [ast.literal_eval(l.split('#', 1)[0].split('=', 1)[1].strip()) for l in open(sys.argv[1]) if l.split('#', 1)[0].split('=', 1)[0].strip() == 'method']; m = v[0] if v else ['snapshot_tree_ilp']; m = [m] if isinstance(m, str) else m; print(m[int(sys.argv[2])])" "$CONFIG" "$SLURM_ARRAY_TASK_ID")
    echo "array task ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}: method ${METHOD}"
    METHOD_ARGS=(--method "$METHOD")
fi
python -u single_cell_inference/snapshot_tree_ilp.py --config "$CONFIG" "${METHOD_ARGS[@]}"
