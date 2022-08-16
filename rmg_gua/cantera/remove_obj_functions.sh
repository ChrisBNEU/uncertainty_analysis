#!/bin/bash
#SBATCH --job-name=parametric_uncertainty
#SBATCH --time=01:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mail-user=$USER@northeastern.edu 
#SBATCH --mail-type=FAIL,END


find /scratch/blais.ch/methanol_results_2022_05_09/ -mindepth 3 -maxdepth 3 -name "*objective_function*" -delete -print
