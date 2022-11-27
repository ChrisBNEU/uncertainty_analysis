import job_manager
import time
import glob
import os
import numpy as np 


from uncertainty_rmg.generate_perturbed_files import generate_perturbed_files
from uncertainty_rmg.make_slurm_scripts import make_slurm_scripts
from uncertainty_cantera.Spinning_basket_reactor.make_slurm_analysis_scripts import make_slurm_analysis_scripts

# specify inputs 
M = 5000 # number of runs of RMG to do
N = 500 # number of slurm jobs to run at once
unc_folder = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
conda_path = unc_folder + "conda/"
RMG_base_folder = unc_folder + "RMG-Py/"
RMG_db_folder = unc_folder + "RMG-database/"
conda_path = unc_folder + "conda/"
expt_yaml_file = ""
# output_path  = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/"
rmg_unc_scripts_folder = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_rmg/"
output_path = "/scratch/blais.ch/methanol_results_2022_10_010_all_params_27_spec/"

# output_csv_name = "ct_analysis.csv"
output_csv_name = "ct_analysis.csv"

###############################################################################
# run all cantera scripts
###############################################################################

# make all of the cantera slurm analysis scripts
working_dir = output_path

make_slurm_analysis_scripts(
    unc_folder, 
    working_dir, 
    conda_path,
    output_csv_name, 
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
