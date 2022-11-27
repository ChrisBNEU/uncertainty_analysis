import job_manager
import time
import glob
import os
import numpy as np 
import yaml

from rmg_gua.rmg.generate_perturbed_files import generate_perturbed_files
from rmg_gua.rmg.make_slurm_scripts import make_slurm_scripts
from rmg_gua.rmg.clean_rmg_db import clean_rmg_db
from rmg_gua.gua_cantera.Spinning_basket_reactor.make_slurm_analysis_scripts import make_slurm_analysis_scripts
from rmg_gua.global_unc.create_uqtk_files import create_files

# specify inputs 
M = 5000 # number of runs of RMG to do
N = 500 # number of slurm jobs to run at once
unc_folder = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
conda_path = os.path.join(unc_folder, "conda")
RMG_base_folder = os.path.join(unc_folder, "RMG-Py")
RMG_db_folder = os.path.join(unc_folder,"RMG-database")
gua_dir = os.path.join(unc_folder, "rmg_gua", "global_unc")
expt_yaml_file = os.path.join(unc_folder, "rmg_gua", "cantera" "all_experiments_reorg.yaml")

rmg_unc_scripts_folder = os.path.join(unc_folder,"rmg_gua","rmg")
output_path = "/scratch/blais.ch/methanol_results_2022_10_06_5000_runs_all_params_27_spec/"
pert_yaml = os.path.join(unc_folder, "rmg_gua", "example", "perturb_groups.yaml")
reference_input = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/rmg/input.py"

# make output folder if it doesn't exist
if not os.path.exists(output_path):
    os.mkdir(output_path)

# define the stuff we want to perturb per yaml file
with open(pert_yaml, 'r') as file:
    perts = yaml.safe_load(file)

##############################################################################
# run all RMG scripts
##############################################################################

# remove files from previous executions
clean_rmg_db()

generate_perturbed_files(
    RMG_db_folder,
    output_path,
    perts,
    M=M,
)

make_slurm_scripts(
    RMG_base_folder, 
    RMG_db_folder,
    unc_folder,
    output_path,
    conda_path,
    rmg_unc_scripts_folder,
    reference_input,
    pert_yaml,
    N=N,
    M=M,
)

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

    
# run uncertainty
create_files(
    output_path,
    gua_dir, 
    M,
)

working_dir = output_path + "rmg_run_scripts/"
print("running GSA scripts")
baseline_job = os.path.join(unc_folder, "rmg_gua", "global_unc", "run_workflow.sh")

# submit baseline: 
gsa_job = job_manager.SlurmJob()
my_cmd = f'sbatch {baseline_job}'
print(my_cmd)
gsa_job.submit(my_cmd)
