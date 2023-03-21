# test the sbr class 
import os
import sys
import yaml
if os.path.exists("/work"):
    prefix = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
    rmg_prefix = prefix
else:
    prefix = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/"
    rmg_prefix = "/Users/blais.ch/Documents/_01_code/RMG_env_1"
sys.path.append(prefix)

from rmg_gua.gua_cantera.Spinning_basket_reactor.sbr import MinSBR
file_path = os.path.join(prefix, "rmg_gua","baseline","cantera","chem_annotated.cti")
kin_par_path = os.path.join(prefix, "rmg_gua","gua_peuqse","02_opt_logp_fam_be","ct_initial_small.yaml")
expt_condts = os.path.join(prefix, "rmg_gua","gua_cantera","experiments_reorg_onlyopt.yaml")

with open(expt_condts, 'r') as file:
    data = yaml.safe_load(file)
conditions = data[2]

test_sbr = MinSBR(
    file_path,
    reac_config=conditions,
    rtol=1.0e-11,
    atol=1.0e-22,
)
ck_matches_old,ck_matches, ck_no_match = test_sbr.test_change_reactions()
len(ck_matches_old)