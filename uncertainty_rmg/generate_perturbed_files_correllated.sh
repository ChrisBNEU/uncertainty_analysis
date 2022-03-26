#!/bin/bash
#SBATCH --job-name=perturb_params
#SBATCH --partition=short
#SBATCH --exclude=c5003
#SBATCH --error=GenErr.log

source activate ../conda
python-jl ./generate_perturbed_files_correllated.py

