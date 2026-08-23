#!/bin/bash
#SBATCH --job-name=data_generation
#SBATCH --nodes=1
#SBATCH --ntasks=1
# networks are generated in parallel, one per cpu (see --n_processes below), so this is the knob that
# sets the generation parallelism. Memory is shared across those processes, each holding one network and
# its matrices, so raise --mem alongside it.
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --mail-type=FAIL,END        # send email when job begins
#SBATCH --mail-user=abronner@princeton.edu

#SBATCH --output=/scratch/gpfs/MONA/abronner/slurm_logs/%j_log.txt
#SBATCH --error=/scratch/gpfs/MONA/abronner/slurm_logs/%j_error.txt

module load anaconda3/2024.6
module load gurobi/8.0.1
conda activate torch_cobra

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"             # must be the project root

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
# Make Gurobi honor the allocation
EXP="${1:?Usage: sbatch run_experiment.sbatch <experiment_name>}" # change per experiment
echo "EXP=$EXP"
echo "EXP=$EXP" >&2

# passed on the command line rather than set per config: configargparse gives the command line priority
# over the config file, so the process count always matches what SLURM actually allocated and no
# experiment config has to be edited when --cpus-per-task changes.
python -m synthetic_data_generation.generate \
    -c configs_for_experiments/data_generation/${EXP}.txt \
    --n_processes "${SLURM_CPUS_PER_TASK:-1}"
