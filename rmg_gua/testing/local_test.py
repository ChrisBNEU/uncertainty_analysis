import yaml
from rmg_gua.gua_rms.sbr import rms_sbr
from rmg_gua.expand_tree.expand_rules import make_new_rule

rerun_sbr = True
if rerun_sbr:
    # testing for local computer before using discovery
    file_path = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/baseline/rms/chem25.rms"
    rmg_db_folder= "/Users/blais.ch/Documents/_01_code/RMG_env_1/RMG-database/"
    expt_condts = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/gua_cantera/all_experiments_reorg_sbr.yaml"

    with open(expt_condts, 'r') as file:
        data = yaml.safe_load(file)

    # pick just one experiment for example 
    conditions = data[8]
    print("run num", conditions['run_num'])
    test_sbr = rms_sbr(
        file_path,
        reac_config = conditions,
        rtol=1.0e-11,
        atol=1.0e-22,
        )
    results = test_sbr.run_simulation()
    


# model_path = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/baseline"
# rmg_db_folder = "/Users/blais.ch/Documents/_01_code/RMG_env_1/RMG-database/"
# make_new_rule(model_path, rmg_db_folder)
