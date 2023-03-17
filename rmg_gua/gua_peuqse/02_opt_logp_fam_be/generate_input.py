import sys
import os 
repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

sys.path.append(repo_dir)
from rmg_gua.gua_cantera.Spinning_basket_reactor.make_peuq_config import *



rmg_path = os.path.dirname(os.environ["RMGPY"])
results_path = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "02_opt_logp_fam_be", "config")
if not os.path.exists(results_path):
    os.mkdir(results_path)


make_rmg_reac_config(
    rmg_path=rmg_path,
    results_path=results_path, 
)   
make_ck_reac_config(results_path=results_path, trim_rules=True)
make_be_peuq_input(results_path=results_path)
make_be_config(results_path=results_path)

