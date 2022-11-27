#!/bin/bash
#SBATCH --job-name=gsa_input
#SBATCH --time=01:00:00
#SBATCH --error=uq_error.log
#SBATCH --output=uq_output.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mail-user=$USER@northeastern.edu 
#SBATCH --mail-type=FAIL,END

# load in the initialization script
source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda
python -u create_uqtk_files.py