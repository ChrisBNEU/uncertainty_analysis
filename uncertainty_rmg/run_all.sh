#!/bin/bash
#SBATCH --job-name=rmg_runs_0-19.sh
#SBATCH --error=run_all_error.log
#SBATCH --output=run_all_output.log
#SBATCH --nodes=1
#SBATCH --partition=west,short
#SBATCH --exclude=c5003
#SBATCH --mem=8Gb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --wait


srun generate_perturbed_files_correllated.py
srun copy_rmg_database_correllated.py
srun run_slurm_scripts_correllated.py


