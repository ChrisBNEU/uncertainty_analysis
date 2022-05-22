#!/bin/bash
#SBATCH --job-name=rmg_base.sh
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --nodes=1
#SBATCH --partition=west
#SBATCH --exclude=c5003
#SBATCH --mem=10Gb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4

export PYTHONPATH="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-Py/:$PYTHONPATH"
export PATH="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-Py/:$PATH"
export RMG="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-Py/rmg.py"
source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda/
python-jl $RMG input.py