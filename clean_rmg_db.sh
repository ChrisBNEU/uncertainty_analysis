#!/bin/bash
#SBATCH --job-name=clean_db
#SBATCH --time=1:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --exclude=c5003
#SBATCH --mem=10Gb
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mail-user=$USER@northeastern.edu 
#SBATCH --mail-type=FAIL,END


source activate ./conda
source ./rmg_gua/set_path_rmg.sh
python -u rmg_gua/rmg/clean_rmg_db.py