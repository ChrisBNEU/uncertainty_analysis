#!/bin/bash
#SBATCH --job-name=uqtk
#SBATCH --time=00:30:00
#SBATCH --error=uq_error.log
#SBATCH --output=uq_output.log
#SBATCH --nodes=1
#SBATCH --partition=west
#SBATCH --exclude=c5003
#SBATCH --mem=5Gb
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mail-user=$USER@northeastern.edu 
#SBATCH --mail-type=FAIL,END

source activate uqtk
source activate_uqtk.sh 
sh workflow_meoh_plots.x
