# clean_rmg_db.py
# cleans out all rules.py, metal.py, 
import glob
import os
import logging
import time

def clean_rmg_db():
    unc_folder = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
    conda_path = os.path.join(unc_folder, "conda")
    RMG_base_folder = os.path.join(unc_folder, "RMG-Py")
    RMG_db_folder = os.path.join(unc_folder,"RMG-database")
    expt_yaml_file = os.path.join(unc_folder, "rmg_gua", "cantera" "all_experiments_reorg.yaml")

    rmg_unc_scripts_folder = os.path.join(unc_folder,"rmg_gua","rmg")
    output_path = "/scratch/blais.ch/methanol_results_2022_08_16_small/"
    pert_yaml = os.path.join(unc_folder, "rmg_gua", "example", "perturb_groups.yaml")
    reference_input = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/rmg/input.py"

    # not very robust but does the job,
    delete_patterns = [
        f'{RMG_db_folder}/input/kinetics/families/*/rules_[0-9][0-9][0-9][0-9].py',    
        f'{RMG_db_folder}/input/kinetics/libraries/Surface/*/*/reactions_[0-9][0-9][0-9][0-9].py',
        f'{RMG_db_folder}/input/thermo/libraries/surfaceThermoPt111_[0-9][0-9][0-9][0-9].py',
        f'{RMG_db_folder}/input/thermo/groups/adsorptionPt111_[0-9][0-9][0-9][0-9].py',
        f'{RMG_db_folder}/input/surface/libraries/metal_[0-9][0-9][0-9][0-9].py',
    ]
    delete_files = []
    for pattern in delete_patterns:
        if len(glob.glob(pattern)) > 0:
            for name in glob.glob(pattern):
                delete_files.append(name)
        else: 
            logging.error(f"could not find any files matching {pattern}")

    # time it and show feeddback, since this part takes a while 
    st = time.time()
    counter = 0
    for name in delete_files:
        if os.path.isfile(name): # this makes the code more robust
            os.remove(name)
            counter +=1
            if counter %1000 == 0:
                elapsed = time.time() - st
                print(f"took {elapsed} sec to del {counter} files")
                print(f"currently deleting {name}")


if __name__ == "__main__":
    clean_rmg_db()