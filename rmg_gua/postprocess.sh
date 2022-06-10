#!/bin/bash
#SBATCH --job-name=postprocess
#SBATCH --time=60:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --nodes=1
#SBATCH --partition=west
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mail-user=$USER@northeastern.edu 
#SBATCH --mail-type=FAIL,END

source activate ./conda
python postprocess_full_model.py