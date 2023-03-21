#!/bin/bash
source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/peuquse_env

# load modules
module load gcc/10.1.0
module load openmpi/4.1.2-gcc10.1

# python-jl run_peuqse.py
export MYMPI="/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/peuquse_env/bin/mpiexec"

echo "MPI is $MYMPI"
# python run_peuqse.py

