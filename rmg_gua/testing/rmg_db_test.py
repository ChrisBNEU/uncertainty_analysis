import unittest
import os
import yaml
import glob
import job_manager
from rmg_gua.rmg.generate_perturbed_files import generate_perturbed_files
from rmg_gua.rmg.make_slurm_scripts import make_slurm_scripts
from rmg_gua.cantera.Spinning_basket_reactor.make_slurm_analysis_scripts import make_slurm_analysis_scripts

# specify inputs 
M = 10 # number of runs of RMG to do
N = 10 # number of slurm jobs to run at once
unc_folder = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
conda_path = os.path.join(unc_folder, "conda")
RMG_base_folder = os.path.join(unc_folder, "RMG-Py")
RMG_db_folder = os.path.join(unc_folder,"RMG-database")
expt_yaml_file = os.path.join(unc_folder, "rmg_gua", "examples" "all_experiments_reorg.yaml")


rmg_unc_scripts_folder = os.path.join(unc_folder,"rmg_gua","rmg")
output_path  = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/"
pert_yaml = os.path.join(unc_folder, "rmg_gua", "example", "perturb_groups.yaml")
reference_input = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/example/input.py"

# define the stuff we want to perturb per yaml file
with open(pert_yaml, 'r') as file:
    perts = yaml.safe_load(file)

##############################################################################
# run all RMG scripts
##############################################################################

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
    N=N,
    M=M,
)

working_dir = output_path + "rmg_run_scripts/"
print("Collecting SLURM scripts")
slurm_scripts = glob.glob(os.path.join(working_dir, "rmg_runs_*.sh"))
baseline_job = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/run_baseline.sh"
slurm_scripts.sort()

# submit baseline: 
rmg_job = job_manager.SlurmJob()
my_cmd = f'sbatch {baseline_job}'
print(my_cmd)
rmg_job.submit(my_cmd)

for i, script in enumerate(slurm_scripts):
    print(f"{i}/{len(slurm_scripts)}\tRunning job {script}")

    rmg_job = job_manager.SlurmJob()
    my_cmd = f'sbatch {script}'
    print(my_cmd)
    rmg_job.submit(my_cmd)

    # wait for job
    # how do we get wait_all to acept a non-array job?
    rmg_job.wait_all()
    print(i)
    


