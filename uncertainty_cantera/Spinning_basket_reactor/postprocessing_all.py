import os
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt

import seaborn
plt.style.use('seaborn-white')

import matplotlib.cm as cm

###############################################################################
# update all objective functions 
###############################################################################

def postprocess_cantera_results(
    output_path
    ):
    
    results = glob.glob(os.path.join(output_path, "run_*/cantera/ct_analysis.csv"))

    for file in results:

        rmg_model_folder = file.replace('cantera/ct_analysis.csv', '')
        cti_file_path = rmg_model_folder = file.replace('ct_analysis.csv', '')
        df = pd.read_csv(file)

        df_graaf = df[(df['T (K)'] < 518) & (df['experiment'] == 'graaf_1988')]
        
        # df['log10(RMG/graaf) MeOH TOF'] = np.log10(max(1e-9, df['RMG MeOH TOF 1/s'] / df_graaf['graaf MeOH TOF 1/s']))
        df_graaf['log10(RMG/graaf) MeOH TOF'] = np.log10((df_graaf['RMG MeOH TOF 1/s'].div(df_graaf['graaf MeOH TOF 1/s']).replace(0, 1e-9)))
        df_graaf['log10(RMG/graaf) H2O TOF'] = np.log10((df_graaf['RMG H2O TOF 1/s'].div(df_graaf['graaf H2O TOF 1/s'])).replace(0, 1e-9))
        df_graaf['log10(RMG/graaf) TOF'] = 0.5 * ( df_graaf['log10(RMG/graaf) MeOH TOF'] + df_graaf['log10(RMG/graaf) H2O TOF'])

        obj_func = df_graaf['log10(RMG/graaf) TOF'].sum()
        print("objective function: ", obj_func)

        # this is naive, but currently saving the objective function to a text file 
        # so we can parse all of them after. 
        obj_func_file = csv_path = os.path.join(rmg_model_folder, "objective_function_log.txt")
        with open(obj_func_file, "w") as f:
            f.write(cti_file_path + ":" + str(obj_func))

output_path  = "/scratch/blais.ch/methanol_results_2022_04_06/"
postprocess_cantera_results(output_path)