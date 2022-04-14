import job_manager
import time
import glob
import os
# import from rmg uncertainty folder
from uncertainty_rmg.generate_perturbed_files import generate_perturbed_files
from uncertainty_rmg.copy_rmg_database import copy_rmg_database
from uncertainty_rmg.make_slurm_scripts import make_slurm_scripts

# import from cantera uncertainty folder
from uncertainty_cantera.Spinning_basket_reactor.make_slurm_analysis_scripts import make_slurm_analysis_scripts

# specify inputs 
M = 5000 # number of runs of RMG to do
N = 100 # number of slurm jobs to run at once
unc_folder = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
conda_path = unc_folder + "conda/"
RMG_base_folder = unc_folder + "RMG-Py/"
RMG_db_folder = unc_folder + "RMG-database/"
conda_path = unc_folder + "conda/"
expt_yaml_file = ""
# output_path  = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/"
output_path = "/scratch/blais.ch/methanol_results_2022_04_06/"

###############################################################################
# run all RMG scripts
###############################################################################

generate_perturbed_files(
    RMG_db_folder,
    unc_folder,
    M=M,
)
copy_rmg_database(
    RMG_db_folder,
    output_path,
    N=N,
)
make_slurm_scripts(
    RMG_base_folder, 
    RMG_db_folder,
    output_path,
    conda_path,
    N=N,
    M=M,
)

# working_dir = os.path.join(os.getcwd(),"uncertainty_output_folder/rmg_run_scripts/")
working_dir = output_path + "rmg_run_scripts/"
print("Collecting SLURM scripts")
slurm_scripts = glob.glob(os.path.join(working_dir, "rmg_runs_*.sh"))

slurm_scripts.sort()

for i, script in enumerate(slurm_scripts):
    print(f"{i}/{len(slurm_scripts)}\tRunning job {script}")

    rmg_job = job_manager.SlurmJob()
    my_cmd = f'sbatch {script}'
    print(my_cmd)
    rmg_job.submit(my_cmd)

    # wait for job
    rmg_job.wait_all()
    print(i)


# ###############################################################################
# # run all cantera scripts
# ###############################################################################

# make all of the cantera slurm analysis scripts
working_dir = output_path

make_slurm_analysis_scripts(
    unc_folder, 
    working_dir, 
    conda_path,
    M=M, 
    N=N,
    )



print("Collecting SLURM scripts")
slurm_scripts = glob.glob(os.path.join(working_dir, "ct_run_scripts/ct_runs_*.sh"))
print(working_dir)
for i, script in enumerate(slurm_scripts):
    print(f"{i}/{len(slurm_scripts)}\tRunning job {script}")

    rmg_job = job_manager.SlurmJob()
    my_cmd = f'sbatch {script}'
    print(my_cmd)
    rmg_job.submit(my_cmd)

    # wait for job
    rmg_job.wait_all()

# parse the objective functions and get the minimum
