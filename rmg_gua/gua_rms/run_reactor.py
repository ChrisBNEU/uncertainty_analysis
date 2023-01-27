import os
import sys
import numpy as np
import pandas as pd
import itertools
import yaml
import cantera as ct
from multiprocessing import Pool
from sbr import rms_sbr
# import time

if len(sys.argv) < 2:
    raise ValueError("Incorrect usage. Must pass the cantera model file as an argument to this analysis script")

if not os.path.exists(sys.argv[1]):
    raise OSError(f"Path to the model file does not exist: {sys.argv[1]}")

# if specified, use a different name for the results file. 
if len(sys.argv) == 3:
    output_file_name = sys.argv[2]
else: 
    output_file_name = "rms_analysis.csv"

rms_file_path = sys.argv[1]
rmg_model_folder = os.path.dirname(rms_file_path)
csv_path = os.path.join(rmg_model_folder, output_file_name)

# generate settings array
settings_yaml = '/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_cantera/all_experiments_reorg_sbr.yaml'
with open(settings_yaml, 'r') as f:
    settings = yaml.safe_load(f)

def run_reactor(condts):

    # initialize reactor
    sbr_ss = rms_sbr(
        rms_file_path,
        reac_config = condts,
        rtol=1.0e-11,
        atol=1.0e-22,
    )

    results = sbr_ss.run_reactor_ss_memory()
    return results


# Too much memory? is that why it's slow?
with Pool() as p:
    result = p.map(run_reactor, settings)

df = pd.DataFrame(result)
df.to_csv(csv_path)

# end = time.time()
# print(f"Completed {len(settings)} processes in {end-start} seconds")

# post process results after pool is finished running
# we will only use runs where intraparticle diffusion limitations
# are not an issue, i.e. T < 518K
df_graaf = df[(df['T (K)'] < 518) & (df['experiment'] == 'graaf_1988') & (df['use_for_opt'] == True)]
obj_func = df_graaf['obj_func'].sum()
print("objective function: ", obj_func)

obj_func_log = df_graaf['log10(RMG/graaf) TOF'].sum()
print("objective function log: ", obj_func_log)


# make objective function title have cantera file name for easier id
ct_run_str = output_file_name.replace(".csv", "")

# this is naive, but currently saving the objective function to a text file 
# so we can parse all of them after.
obj_func_file = os.path.join(rmg_model_folder, f"objective_function_{ct_run_str}.txt")
with open(obj_func_file, "w") as f:
    f.write(rms_file_path + ":" + str(obj_func))

# make log_obj_func file
obj_func_file_log = os.path.join(rmg_model_folder, f"objective_function_log_{ct_run_str}.txt")
with open(obj_func_file_log, "w") as f:
    f.write(rms_file_path + ":" + str(obj_func_log))