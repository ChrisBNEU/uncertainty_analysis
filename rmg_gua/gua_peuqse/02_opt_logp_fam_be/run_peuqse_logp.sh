#!/bin/bash
#SBATCH --job-name=peuq_logp
#SBATCH --time=8:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --tasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --mem=200Gb
#SBATCH --ntasks=128
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_peuqse/02_opt_logp_fam_be
source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/peuquse_env

# load modules
module load gcc/10.1.0
module load openmpi/4.1.2-gcc10.1

# python-jl run_peuqse.py
# /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/peuquse_env/bin/mpiexec -n 28 python run_peuqse.py
mpiexec -n 28 -v python run_peuqse.py
# python run_peuqse.py

