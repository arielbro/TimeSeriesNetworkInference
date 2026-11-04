#!/bin/bash
# Submits run_single_cell_tree_ilp.sbatch as a SLURM array with one task per method in the config's method list.
# Run from the project root; needs only a python on the path (the config is parsed without imports). Extra arguments go to sbatch (e.g. --mem=156G).
# Usage: bash configs_for_experiments/slurm_scripts/submit_single_cell_tree_ilp.sh [<config name>] [sbatch options]
set -euo pipefail
EXP="${1:-single_cell_tree_ilp}"
shift || true
N=$(python -c "import ast, sys; v = [ast.literal_eval(l.split('#', 1)[0].split('=', 1)[1].strip()) for l in open(sys.argv[1]) if l.split('#', 1)[0].split('=', 1)[0].strip() == 'method']; m = v[0] if v else ['snapshot_tree_ilp']; m = [m] if isinstance(m, str) else m; print(len(m))" "configs_for_experiments/inference/${EXP}.txt")
sbatch --array=0-$((N - 1)) "$@" configs_for_experiments/slurm_scripts/run_single_cell_tree_ilp.sbatch "$EXP"
