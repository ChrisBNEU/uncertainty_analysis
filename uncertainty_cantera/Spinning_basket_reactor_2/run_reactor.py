import os
import sys
import numpy as np
import pandas as pd
import itertools
import yaml
from multiprocessing import Pool
from sbr import MinSBR
# import time

# if len(sys.argv) < 2:
#     raise ValueError("Incorrect usage. Must pass the cantera model file as an argument to this analysis script")

# if not os.path.exists(sys.argv[1]):
#     raise OSError(f"Path to the cantera model file does not exist: {sys.argv[1]}")


# cti_file_path = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/run_0001/cantera/chem_annotated.cti"
cti_file_path = "/work/westgroup/ChrisB/_01_MeOH_repos/meOH-synthesis/base/cantera/chem_annotated.cti"
# cti_file_path = sys.argv[1]
rmg_model_folder = os.path.dirname(cti_file_path)
csv_path = os.path.join(rmg_model_folder, "ct_analysis.csv")


# generate settings array
settings_yaml = '../all_experiments_reorg_sbr.yaml'
with open(settings_yaml, 'r') as f:
    settings = yaml.safe_load(f)

def run_reactor(condts):

    # initialize reactor
    sbr_ss = MinSBR(
        cti_file_path,
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
df_graaf = df[(df['T (K)'] < 518) & (df['experiment'] == 'graaf_1988')]
obj_func = df_graaf['obj_func'].sum()
print("objective function: ", obj_func)

# this is naive, but currently saving the objective function to a text file 
# so we can parse all of them after. 
obj_func_file = csv_path = os.path.join(rmg_model_folder, "objective_function.txt")
with open(obj_func_file, "w") as f:
    f.write(cti_file_path + ":" + str(obj_func))
