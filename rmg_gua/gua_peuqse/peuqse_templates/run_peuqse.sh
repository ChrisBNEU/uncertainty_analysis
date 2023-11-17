#!/bin/bash
#SBATCH --job-name=par_rem
#SBATCH --time=24:00:00
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --tasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --mem=50Gb
#SBATCH --ntasks=128
#SBATCH --mail-user=blais.ch@northeastern.edu
#SBATCH --mail-type=FAIL,END
module purge
source /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/peuquse_env/bin/activate
module load gcc/10.1.0
module load openmpi/4.1.2-gcc10.1
mpiexec -n 128  -v python run_peuqse.py
