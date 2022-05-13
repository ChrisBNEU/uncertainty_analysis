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
output_path = "/scratch/blais.ch/methanol_results_2022_05_09/"

# output_csv_name = "ct_analysis.csv"
output_csv_name = "ct_analysis_artol.csv"

# define the stuff we want to perturb
# Define the number of perturbations to run
# testing currently, use 10
thermo_libraries = [
    'surfaceThermoPt111',
]

# groups as in the group files (adsorptionpt111, ni111, etc.)
thermo_groups_to_perturb = [
    'adsorptionPt111',
]

# pick which entries to perturb in the kinetics library
# WARNING: does not handle overlap of entries in different libraries
# if 'all', perturb all entries
lib_entries_to_perturb = [
    'all'
]
# pick which kinetics libraries to perturb
kinetics_libraries = [
    'CPOX_Pt/Deutschmann2006_adjusted'
]


# Load the databases
kinetics_families = [  # list the families to perturb
        "Surface_Abstraction",
        "Surface_Abstraction_Beta",
        "Surface_Abstraction_Beta_double_vdW",
        "Surface_Abstraction_Beta_vdW",
        "Surface_Abstraction_Single_vdW",
        "Surface_Abstraction_vdW",
        "Surface_Addition_Single_vdW",
        "Surface_Adsorption_Abstraction_vdW",
        "Surface_Adsorption_Bidentate",
        "Surface_Adsorption_Dissociative",
        "Surface_Adsorption_Dissociative_Double",
        "Surface_Adsorption_Double",
        "Surface_Adsorption_Single",
        "Surface_Adsorption_vdW",
        "Surface_Bidentate_Dissociation",
        "Surface_Dissociation",
        "Surface_Dissociation_Beta",
        "Surface_Dissociation_Beta_vdW",
        "Surface_Dissociation_Double",
        "Surface_Dissociation_Double_vdW",
        "Surface_Dissociation_vdW",
        "Surface_DoubleBond_to_Bidentate",
        "Surface_Dual_Adsorption_vdW",
        "Surface_EleyRideal_Addition_Multiple_Bond",
        "Surface_Migration",
        "Surface_vdW_to_Bidentate",
]

# pass in values as a dictionary to make it easier
perturb_dict = {
    "thermo_libraries" : thermo_libraries,
    "thermo_groups" : thermo_groups_to_perturb,
    "lib_entries_to_perturb" : lib_entries_to_perturb,
    "kinetic_libraries" : kinetics_libraries,
    "kinetics_families" : kinetics_families,
}

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
