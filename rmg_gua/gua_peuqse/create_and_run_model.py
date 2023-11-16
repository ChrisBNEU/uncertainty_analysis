import shutil
import os
import sys
import argparse
import subprocess
repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
print("repo_dir: ", repo_dir)
sys.path.append(repo_dir)


# simple script for creating a new run folder and running an rmg model
def create_run_folder(full_path, parallel=False, local_run=False, overwrite=False):


    rmg_model_folder = os.path.join(full_path, "rmg_model")

    if not os.path.exists(full_path):
        os.makedirs(full_path)
        os.mkdir(rmg_model_folder)

    # copy the rmg model input file
    rmg_input_file = os.path.join(repo_dir, "rmg_gua", "baseline", "input.py")
    shutil.copy(rmg_input_file, rmg_model_folder)

    # copy the run script
    run_script = os.path.join(repo_dir, "rmg_gua", "baseline", "run_rmg.sh")
    run_script_dest = os.path.join(rmg_model_folder, "run_rmg.sh")

    with open(run_script, "r") as f:
        lines = f.readlines()
        with open(run_script_dest, "w") as f2: 
            for line in lines:
                if line.startswith("cd MODEL_PATH"): 
                    line = line.replace("MODEL_PATH", rmg_model_folder)    
                f2.write(line) 


    return rmg_model_folder


def run_rmg_model(model_path):
    # run the model
    os.chdir(model_path)
    subprocess.run(["bash", "run_rmg.sh"])
    os.chdir(os.path.abspath(__file__))
    


if __name__ == "__main__":

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

    #cmd line args
    args = parser.parse_args()
    run_folder = args.folder
    parallel = args.parallel
    local_run = args.local

    # optional arg for scratch path (used on nersc)
    file_path = args.directory
    full_path = os.path.join(file_path, "peuqse_runs", run_folder)

    # make the run folder 
    rmg_model_folder = create_run_folder(full_path, parallel=parallel, local_run=local_run)
    run_rmg_model(rmg_model_folder)