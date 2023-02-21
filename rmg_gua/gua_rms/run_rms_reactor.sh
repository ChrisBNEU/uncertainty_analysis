#!/bin/bash
#SBATCH --job-name=sbr_runs
#SBATCH --time=1:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

source activate ../../conda
source ../set_path_rmg.sh
python-jl run_reactor.py /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/rms/chem53.rms

