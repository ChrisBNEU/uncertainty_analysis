#!/bin/bash
#SBATCH --job-name=copy_rmg_database
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c5003
#SBATCH --mem=8Gb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --error=dberror.log


# Copy the files from the full database to the mostly symbolic one
python ./copy_rmg_database_correllated.py



