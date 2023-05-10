#!/bin/bash
#SBATCH --job-name=makefolder
#SBATCH --time=1:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=10Gb
#SBATCH --ntasks=1
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_peuqse/
source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda/

source /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/set_path_rmg.sh

python-jl setup_run_folder.py 04_logp_opt_linear_checkpoint 0 0
