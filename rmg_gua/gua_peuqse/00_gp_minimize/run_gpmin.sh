#!/bin/bash
#SBATCH --job-name=meohgpmin
#SBATCH --time=5:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=short
#SBATCH --mem=100Gb
#SBATCH --ntasks=1
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_peuqse/00_gp_minimize
source activate /work/westgroup/ChrisB/_06_chevron/03_kme/chevron

# python-jl run_peuqse.py
python -u gp_min_script.py

