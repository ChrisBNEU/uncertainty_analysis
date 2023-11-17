#!/bin/bash
DEST="/scratch/blais.ch/methanol_unc_data/"
FOLDER="peuqse_methanol_runs"

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_peuqse/
source /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda/bin/activate

source /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/set_path_rmg.sh

# first, run rmg model using current db and py
python-jl create_and_run_model.py -p -f $FOLDER -d $DEST

# then, fill with the prerequisite files
python-jl setup_run_folder.py -p -f $FOLDER -d $DEST
