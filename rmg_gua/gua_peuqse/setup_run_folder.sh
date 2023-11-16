#!/bin/bash
#SBATCH --job-name=makefolder
#SBATCH --time=1:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=10Gb
#SBATCH --ntasks=1
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

DEST="/scratch/blais.ch/methanol_unc_data/"
FOLDER="peuqse_methanol_runs"

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_peuqse/
source /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda/bin/activate

source /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/set_path_rmg.sh

# first, run rmg model using current db and py
python-jl create_and_run_model.py -p -f $FOLDER -d $DEST

# then, fill with the prerequisite files
python-jl setup_run_folder.py -p -f $FOLDER -d $DEST