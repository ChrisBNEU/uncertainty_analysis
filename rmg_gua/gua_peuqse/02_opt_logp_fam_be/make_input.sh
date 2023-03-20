#!/bin/bash
#SBATCH --job-name=peuq_logp
#SBATCH --time=1:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=50Gb
#SBATCH --ntasks=1
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_peuqse/02_opt_logp_fam_be
source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda

# set env variables
source /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/set_path_rmg.sh

# generate peuqse input files
python-jl generate_input.py
