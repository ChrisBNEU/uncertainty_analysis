#!/bin/bash
#SBATCH --job-name=peuq_logp
#SBATCH --time=8:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=short
#SBATCH --mem=100Gb
#SBATCH --ntasks=1
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_peuqse/02_opt_logp_fam_be
source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/peuquse_env

# load modules
module load gcc/10.1.0
# module load openmpi/3.1.2

echo "running script"
# python-jl run_peuqse.py
python -u run_peuqse.py

