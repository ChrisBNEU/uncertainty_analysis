#!/bin/bash
#SBATCH --job-name=peuq_run
#SBATCH --time=1:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --tasks-per-node=15
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --mem=50Gb
#SBATCH --ntasks=30
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_peuqse

source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/peuquse_env

# load modules
module load gcc/10.1.0
module load openmpi/4.1.2-gcc10.1

# python-jl run_peuqse.py
/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/peuquse_env/bin/mpiexec -n 28 python run_peuqse.py

