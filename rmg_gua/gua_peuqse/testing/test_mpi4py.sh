#!/bin/bash
#SBATCH --job-name=debug_mpi
#SBATCH --time=00:00:30
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH -c 1
##SBATCH --constraint=ib
#SBATCH --exclusive
#SBATCH --partition=debug
#SBATCH --mem-per-cpu=1Gb
#SBATCH -n 5
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

export SLURM_EXPORT_ENV=All

cd /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_peuqse

source activate ../../conda
source ../set_path_rmg.sh

# load modules
module load gcc/10.1.0
# module load openmpi/4.1.2-gcc10.1
module load openmpi/3.1.2

/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/conda/bin/mpiexec -n 5 python test_mpi4py.py
