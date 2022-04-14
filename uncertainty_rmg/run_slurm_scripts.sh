#!/bin/bash
#SBATCH --job-name=parametric_uncertainty
#SBATCH --time=23:59:59
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mail-user=$USER@northeastern.edu 
#SBATCH --mail-type=FAIL,END

source activate ../conda
python ./run_slurm_scripts_correllated.py

