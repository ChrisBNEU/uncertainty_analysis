#!/bin/bash
#SBATCH --job-name=uqtk
#SBATCH --time=1:00:00
#SBATCH --error=uq_error.log
#SBATCH --output=uq_output.log
#SBATCH --nodes=1
#SBATCH --mem=100Gb
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --mail-user=$USER@northeastern.edu 
#SBATCH --mail-type=FAIL,END

cd ./rmg_gua/global_unc
source activate uqtk
source activate_uqtk.sh 
sh workflow_meoh.x
