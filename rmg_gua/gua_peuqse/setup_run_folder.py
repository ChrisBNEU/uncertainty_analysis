import shutil
import os
import sys
repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
print("repo_dir: ", repo_dir)
sys.path.append(repo_dir)
from rmg_gua.gua_cantera.Spinning_basket_reactor.make_peuq_config import make_rmg_reac_config, make_be_peuq_input, \
    trim_rule_file, make_ck_reac_config, make_be_config
from rmg_gua.gua_peuqse.peuqse_utilities import make_ct_expt_file, make_lookup_dict

# first command line arg is the run folder location 
# second is bool. if true do simple sum squares residual opt
# third is bool. if true do by species perturbations. false, do binding energy perts.
def make_bash_script(full_path, parallel=False, local_run=False):
    """ 
    make the bash/slurm script depending on the options selected
    """
    if parallel:
        run_type = "par"
        # define nodes and total number of tasks you want to run in parallel
        nodes = 2
        tot_tasks = 128
        tpn = int(tot_tasks / nodes)

        # mpi runs need one node extra to control, so # of sub tasks is n-1
        mpi_tasks = tpn - 1

    else: 
        run_type = "ser"
        tpn = 1
        tot_tasks = 1
    
    if local_run: 
        run_loc = "loc"
    else: 
        run_loc = "rem"
    
    bash_path = os.path.join(full_path, "run_peuqse.sh")
    
    output = []
    output.append("#!/bin/bash")
    if not local_run:
        output.append(f"#SBATCH --job-name={run_type}_{run_loc}")
        output.append("#SBATCH --time=24:00:00")
        output.append("#SBATCH --error=error.log")
        output.append("#SBATCH --output=output.log")
        output.append(f"#SBATCH --tasks-per-node={tpn}")
        output.append("#SBATCH --cpus-per-task=1")
        output.append("#SBATCH --partition=short")
        output.append("#SBATCH --mem=50Gb")
        output.append(f"#SBATCH --ntasks={tot_tasks}")
        output.append("#SBATCH --mail-user=blais.ch@northeastern.edu")
        output.append("#SBATCH --mail-type=FAIL,END")

        output.append("source activate /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/peuquse_env")
        output.append("module load gcc/10.1.0")
        output.append("module load openmpi/4.1.2-gcc10.1")
        output.append("rm /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/cantera/chem_annotated.yaml")
        output.append("cti2yaml /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/cantera/chem_annotated.cti")
    else: 
        output.append("source activate peuqse_env")

    
    if parallel: 
        output.append(f"mpiexec -n {tot_tasks} python run_peuqse.py")
    else:
        output.append("python run_peuqse.py")

    with open(bash_path, "w") as f:
        for line in output: 
            f.write(line+"\n")

def make_run_file(full_path, parallel=False):
    """
    make the python file for running peuqse
    """
    run_path = os.path.join(full_path, "run_peuqse.py")
    output = []
    output.append("import copy")
    output.append("import pickle")
    output.append("import shutil")
    output.append("import os")
    output.append("import math")
    output.append("import yaml")
    output.append("import numpy as np")
    output.append("import PEUQSE as PEUQSE")
    output.append("import PEUQSE.UserInput as UserInput")
    output.append("import sys")
    output.append("repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))")
    output.append("sys.path.insert(0, repo_dir)")
    output.append("import rmg_gua.gua_peuqse.ct_simulation as ct_simulation")
    output.append("from rmg_gua.gua_peuqse.setup_peuqse import setup_userinput")
    output.append("project_path = os.path.dirname(os.path.abspath(__file__))")
    output.append("")
    output.append("if __name__ == \"__main__\":")
    output.append("    print(\"running job\")")
    output.append("    # setup our ct_simulation function")
    output.append("    ct_simulation.sim_init(project_path)")
    output.append("    UserInput = setup_userinput(project_path)")
    output.append("    UserInput.model['exportResponses'] = True")
    output.append("    UserInput.parameter_estimation_settings['mcmc_threshold_filter_samples'] = True")
    output.append(f"    UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] = {parallel}")
    output.append("    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0")
    output.append("    UserInput.parameter_estimation_settings['multistart_searchType'] = 'doOptimizeLogP'")
    output.append("    UserInput.parameter_estimation_settings['multistart_passThroughArgs'] = {'method':'Nelder-Mead'}")
    output.append("    UserInput.parameter_estimation_settings['multistart_initialPointsDistributionType'] = 'sobol'")
    output.append("    UserInput.parameter_estimation_settings['multistart_exportLog'] = True")
    output.append("    UserInput.parameter_estimation_settings['multistart_checkPointFrequency'] = 100")
    output.append("    UserInput.parameter_estimation_settings['verbose'] = True")
    output.append("    PE_object = PEUQSE.parameter_estimation(UserInput)")
    output.append("    PE_object.doMultiStart()")
    output.append("    PE_object.createAllPlots()")

    with open(run_path, "w") as f:
        for line in output: 
            f.write(line+"\n")

def make_run_folder(full_path, parallel=False, local_run=False):
    """
    makes a folder for peuqse to run in, with all required inputs
    """
    if not os.path.exists(full_path):
        os.mkdir(full_path)
    else: 
        print("folder already exists, updating files")

    # if we are consistent with env variables this should work
    rmg_path = os.path.dirname(os.environ["RMGPY"])

    # base_path
    base_path = os.path.join(repo_dir, "rmg_gua", "baseline")

    # make the input files into a "config" folder
    config_path = os.path.join(full_path, "config")
    if not os.path.exists(config_path):
        os.mkdir(config_path)
    else:
        print("config folder already exists, updating files")

    # specify run path
    results_path = os.path.join(full_path, "config")

    # make yamls for our expt response and uncertainty
    make_ct_expt_file(results_path=results_path, use_peuq_expts=True)

    # make all input files for reactions
    make_rmg_reac_config(rmg_path=rmg_path, results_path=results_path)
    make_ck_reac_config(results_path=results_path, trim_rules=True)

    # make the value, uncertainty, upper bound, and lower bound files for our parameters
    make_be_peuq_input(results_path = results_path)
    make_be_config(results_path=results_path, return_test_spec=True)

    # make the scripts for running the ct simulation
    make_bash_script(full_path, parallel=parallel, local_run=local_run)
    make_run_file(full_path, parallel=parallel)

    # make the lookup dict
    make_lookup_dict(base_path, results_path=results_path)

if __name__ == "__main__":

    run_folder = str(sys.argv[1])
    parallel = bool(int(sys.argv[2]))
    local_run = bool(int(sys.argv[3]))
    file_path = os.path.dirname(__file__)
    full_path = os.path.join(file_path, "peuqse_runs", run_folder)

    make_run_folder(full_path, parallel=parallel, local_run=local_run)

    


    

        