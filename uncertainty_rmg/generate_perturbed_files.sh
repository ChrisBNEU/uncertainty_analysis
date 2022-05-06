#!/bin/bash
#SBATCH --job-name=perturb_params
#SBATCH --partition=express
#SBATCH --error=GenErr.log

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_rmg
source activate ../conda
python-jl ./generate_perturbed_files.py

