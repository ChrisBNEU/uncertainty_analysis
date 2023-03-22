import sys
import os

from multiprocessing import Pool
import subprocess
import re
import numpy as np
import pandas as pd
from skopt import gp_minimize
from skopt import callbacks
from skopt.callbacks import CheckpointSaver
from skopt import dump, load
import yaml

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(repo_dir)
results_path = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "00_gp_minimize", "config")
from rmg_gua.gua_peuqse.runtime_utilities import get_all_param_lists, make_exp_data_lists

from ct_simulation import simulationFunction
def main():
    # optimize
    # load checkpoint if it exists
    use_skopt_checkpoint = True
    ckpt_path = os.path.join(results_path, "checkpoint.pkl")
    
    # get the peuqse parameters
    param_dict = get_all_param_lists(results_path=results_path)     

    # load the exp data yaml
    expt_data_path = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "02_opt_logp_fam_be", "expt_data.yaml")
    with open(expt_data_path, "r") as f:
        expt_data = yaml.load(f, Loader=yaml.FullLoader)
        
    # load uncertainty yaml. exyperimenting with different layout for X. 
    expt_unc_path = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "02_opt_logp_fam_be", "expt_unc.yaml")
    with open(expt_unc_path, "r") as f:
        expt_unc = yaml.load(f, Loader=yaml.FullLoader)
    
    data_path = os.path.join(repo_dir, "rmg_gua", "gua_cantera")
    _, x_data, y_data, y_unc = make_exp_data_lists(data_path, results_path, use_count=3)
    
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    y_unc = np.array(y_data)
    print(f"length is {len(x_data)} in main")
                             
    
    if os.path.exists(ckpt_path):
        res = load(ckpt_path)
        x0 = res.x_iters
        y0 = res.func_vals
    else:
        x0 = param_dict["value_list"]
        y0 = None
    
    # make list of tuples for the ranges
    param_range = []
    for index, value in enumerate(param_dict["value_list"]):
        assert value <= param_dict["upper_list"][index]
        assert value >= param_dict["lower_list"][index]
        ub = float(param_dict["upper_list"][index])
        lb = float(param_dict["lower_list"][index])
        param_range.append((lb, ub))

    checkpoint_saver = CheckpointSaver(ckpt_path, compress=9)
    res = gp_minimize(
        simulationFunction,
        # we initially guessed to match up points that starting cov was 0.6. use this as upper limit
        param_range,
        n_calls=5000,
        x0=x0,
        y0=y0,
        callback=[checkpoint_saver],
        # n_initial_points = 64,
        # initial_point_generator="sobol"
    )
    dump(res, 'result.pkl')
    print("x*=%.2f f(x*)=%.2f" % (res.x[0], res.fun))

        
    print("model complete")

if __name__ == "__main__":
    main()