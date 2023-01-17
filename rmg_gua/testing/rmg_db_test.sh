#!/bin/bash
#SBATCH --job-name=parametric_uncertainty
#SBATCH --time=120:00:00
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

cd ../../
source activate ./conda
python rmg_gua/testing/rmg_db_test.py