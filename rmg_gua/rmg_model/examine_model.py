import os
import sys
import numpy as np
import pandas as pd
import itertools
import yaml
import multiprocessing
# from multiprocessing import Pool

repo_dir = os.path.dirname(os.path.dirname(os.path.abspath("")))
sys.path.append(repo_dir)
import time 
impt1 = time.time()
from rmg_gua.gua_rms.sbr import rms_sbr
impt2 = time.time()
print(f"import of julia packages took {impt2-impt1} seconds")

# if specified, use a different name for the results file. 
if len(sys.argv) == 3:
    output_file_name = sys.argv[2]
else: 
    output_file_name = "rms_analysis.csv"

rms_file_path = sys.argv[1]
rmg_model_folder = os.path.dirname(rms_file_path)
csv_path = os.path.join(rmg_model_folder, output_file_name)

# generate settings array

settings_yaml =os.path.join(repo_dir, "rmg_gua", "gua_cantera", "experiments_reorg_onlyopt.yaml")
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

    results = sbr_ss.run_simulation()
    return results


# Too much memory? is that why it's slow?
# multiprocessing.set_start_method("spawn")
# with multiprocessing.Pool() as p:
#     result = p.map(run_reactor, settings)

result = []
# workaround: just run in serial
count = 0
ttot1 = time.time()
for condition in settings: 
    print("running ", count)
    
    t1 = time.time()
    res = run_reactor(condition)
    result.append(res)
    t2 = time.time()
    
    print(f"process took {t2-t1} seconds")

ttot2 = time.time()
print(f"all reactors took {ttot2-ttot1} seconds")

df = pd.DataFrame(result)
df.to_csv(csv_path)

# end = time.time()
# print(f"Completed {len(settings)} processes in {end-start} seconds")

# post process results after pool is finished running
# we will only use runs where intraparticle diffusion limitations
# are not an issue, i.e. T < 518K
df_graaf = df[(df['T (K)'] < 518) & (df['experiment'] == 'graaf_1988') & (df['use_for_opt'] == True)]

# obj_func = df_graaf['obj_func'].sum()
obj_func = df_graaf['Sum Error %'].sum()
print("objective function: ", obj_func)

# make objective function title have cantera file name for easier id
run_str = output_file_name.replace(".csv", "")

# this is naive, but currently saving the objective function to a text file 
# so we can parse all of them after.
obj_func_file = os.path.join(rmg_model_folder, f"objective_function_{run_str}.txt")
with open(obj_func_file, "w") as f:
    f.write(rms_file_path + ":" + str(obj_func))