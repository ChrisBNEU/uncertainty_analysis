import job_manager
import time
import glob
import os

# specify inputs 
M = 10 # number of runs of RMG to do
N = 50 # number of slurm jobs to run at once
unc_folder = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/
RMG_base_folder = unc_folder + "RMG-Py/"
RMG_db_folder = unc_folder + "RMG-database/"
conda_path = unc_folder + "conda/"
expt_yaml_file = ""
output_path  = unc_folder + "uncertainty_output_folder/"

###############################################################################
# run all RMG scripts
###############################################################################

generate_perturbed_files_correllated()
copy_rmg_database_correllated()

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


###############################################################################
# run all cantera scripts
###############################################################################

# make all of the cantera slurm analysis scripts
make_slurm_analysis_scripts(unc_folder, working_dir, M=10, N=50):

working_dir = output_path + "ct_run_scripts/"

print("Collecting SLURM scripts")
slurm_scripts = glob.glob(os.path.join(working_dir, "ct_runs_*.sh"))

for i, script in enumerate(slurm_scripts):
    print(f"{i}/{len(slurm_scripts)}\tRunning job {script}")

    rmg_job = job_manager.SlurmJob()
    my_cmd = f'sbatch {script}'
    print(my_cmd)
    rmg_job.submit(my_cmd)

    # wait for job
    rmg_job.wait_all()