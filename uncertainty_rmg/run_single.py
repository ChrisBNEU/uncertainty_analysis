# from make_slurm_scripts 
from uncertainty_rmg.make_slurm_scripts import make_slurm_scripts
from uncertainty_rmg.copy_rmg_database import copy_rmg_database

def run_single(
    RMG_base_folder, 
    RMG_db_folder,
    output_path,
    conda_path,
    rmg_unc_scripts_folder,
    perturb_dict,
    N=10, #number of runs to do at once
    M=10, # number of runs to do total
    ):

    make_slurm_scripts(
        RMG_base_folder, 
        RMG_db_folder,
        output_path,
        conda_path,
        rmg_unc_scripts_folder,
        N=10, #number of runs to do at once
        M=10, # number of runs to do total
    )
    copy_rmg_database(
        RMG_db_folder, 
        output_path,
        perturb_dict,
        N,
    )


