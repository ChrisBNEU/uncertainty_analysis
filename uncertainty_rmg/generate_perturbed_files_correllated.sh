#!/bin/bash
#SBATCH --job-name=perturb_params
#SBATCH --partition=short
#SBATCH --exclude=c5003

source activate ../conda
python-jl /scratch/westgroup/methanol/meOH-synthesis/perturbed_runs/generate_perturbed_files_correllated.py

