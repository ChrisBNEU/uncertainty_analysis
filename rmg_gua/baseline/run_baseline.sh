#!/bin/bash
#SBATCH --job-name=rmg_base.sh
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --exclude=c5003
#SBATCH --mem=10Gb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4

cd "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/"
export PYTHONPATH="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-Py/:$PYTHONPATH"
export PATH="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-Py/:$PATH"
export RMG="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-Py/rmg.py"
source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda/

# run RMG job
# python-jl $RMG input.py

MECH_FILE=$(ls "./rms/" | tail -1)
WRKDIR=$PWD

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_rms/
# run cantera
python-jl /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_rms/run_reactor.py $WRKDIR/rms/$MECH_FILE rms_analysis.csv
              

