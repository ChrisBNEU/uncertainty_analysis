import shutil
import os
import sys
import argparse
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

        output.append("module purge")
        output.append("source /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/peuquse_env/bin/activate")
        output.append("module load gcc/10.1.0")
        output.append("module load openmpi/4.1.2-gcc10.1")
        # output.append(f"rm {os.path.join(full_path,'cantera/chem_annotated.yaml')}")
        # output.append("cti2yaml /work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/s/cantera/chem_annotated.cti")
    else: 
        output.append("source activate peuqse_env")

    
    if parallel: 
        output.append(f"mpiexec -n {tot_tasks}  -v python run_peuqse.py")
    else:
        output.append("python run_peuqse.py")

    with open(bash_path, "w") as f:
        for line in output: 
            f.write(line+"\n")

def make_run_file(full_path, parallel=False, run_type="opt"):
    """
    make the python file for running peuqse
    """
    if run_type == "ssr":
        filename = "run_peuqse_optlogp.py"
    elif run_type == "mcmc":
        filename = "run_peuqse_mcmc.py"
    else:
        raise ValueError("run_type must be opt or mcmc")
    
    template_path = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "peuqse_templates", filename)
    run_path = os.path.join(full_path, "run_peuqse.py")

    template_path_slurm = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "peuqse_templates", "run_peuqse.sh")
    run_path_slurm = os.path.join(full_path, "run_peuqse.sh")

    # copy slurm script
    shutil.copy(template_path_slurm, run_path_slurm)

    # copy peuqse script
    with open(template_path, "r") as f:
        lines = f.readlines()
        with open(run_path, "w") as f2:
            for line in lines: 
                if "__REPO_DIR__" in line:
                    line = line.replace("__REPO_DIR__", f"'{repo_dir}'")
                if "__PARALLEL__" in line: 
                    line = line.replace("__PARALLEL__", str(parallel))

                f2.write(line)

def make_run_folder(full_path, parallel=False, local_run=False, run_type="opt"):
    """
    makes a folder for peuqse to run in, with all required inputs
    """
    if not os.path.exists(full_path):
        os.makedirs(full_path)
    else: 
        print("folder already exists, updating files")

    # if we are consistent with env variables this should work
    rmg_path = os.path.dirname(os.environ["RMGPY"])

    # base_path
    # base_path = os.path.join(repo_dir, "rmg_gua", "baseline")
    base_path = os.path.join(full_path, "rmg_model")

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
    make_ck_reac_config(results_path=results_path, trim_rules=True, model_path=base_path)

    # make the value, uncertainty, upper bound, and lower bound files for our parameters
    make_be_peuq_input(results_path = results_path)
    make_be_config(results_path=results_path, return_test_spec=True, model_path=base_path)

    # make the scripts for running the ct simulation
    make_bash_script(full_path, parallel=parallel, local_run=local_run)
    make_run_file(full_path, parallel=parallel, run_type=run_type)

    # make the lookup dict
    make_lookup_dict(base_path, results_path=results_path)

if __name__ == "__main__":

    # run_folder = str(sys.argv[1])
    # parallel = bool(int(sys.argv[2]))
    # local_run = bool(int(sys.argv[3]))

    # if len(sys.argv) > 4:
    #     file_path = str(sys.argv[4])
    # else: 
    #     file_path = os.path.dirname(__file__)
    
    # full_path = os.path.join(file_path, "peuqse_runs", run_folder)

    # make_run_folder(full_path, parallel=parallel, local_run=local_run)


     # cmd line args
    parser = argparse.ArgumentParser()

    # current folder
    curr_folder = os.path.dirname(os.path.abspath(__file__))

    # -f run_folder -o opt_type -p param_type -d directory -n nersc
    # e.g. -f 01_test_folder -p -d /pscratch/sd/c/chrisjb/cpox_peuqse_runs
    parser.add_argument(
        "-f", "--folder", type=str, help="Run folder name"
        )

    parser.add_argument(
        "-r", "--runtype", type=str, help="Optimization type (opt or mcmc)", 
        default='ssr'
        )

    parser.add_argument(
        "-p", "--parallel", 
        help="run in parallel or serial. default is serial",
        default=False, action='store_true'
        )

    parser.add_argument(
        "-l", "--local", 
        help="run locally or on hpc. default is hpc run",
        default=False, action='store_true'
        )
    

    parser.add_argument(
        "-d", "--directory", type=str, help="directory to make run folder in", 
        default=curr_folder
        )

    # parser.add_argument(
    #     "-n", "--nersc", help="Run on nersc, y/n", default=False, action='store_true'
    #     )

    #cmd line args
    args = parser.parse_args()
    run_folder = args.folder
    parallel = args.parallel
    local_run = args.local
    run_type = args.runtype

    # optional arg for scratch path (used on nersc)
    file_path = args.directory
    full_path = os.path.join(file_path, "peuqse_runs", run_folder)
    
    # optional arg to use nersc scripts
    # nersc = args.nersc

    # make the run folder 
    make_run_folder(full_path, parallel=parallel, local_run=local_run, run_type=run_type)

    


    

        