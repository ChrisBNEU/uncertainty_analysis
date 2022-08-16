#!/bin/bash
#SBATCH --job-name=rmg_base.sh
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --nodes=1
#SBATCH --partition=west
#SBATCH --exclude=c5003
#SBATCH --mem=10Gb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4

# activate conda environment
source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda/
python /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_cantera/Spinning_basket_reactor/run_reactor.py /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_baseline/cantera/chem_annotated.cti ct_analysis.csv
              