#!/usr/bin/env python
# coding: utf-8

import glob
import pandas as pd
import numpy as np
import os
import re
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

unc_folder = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
conda_path = unc_folder + "conda/"
RMG_base_folder = unc_folder + "RMG-Py/"
RMG_db_folder = unc_folder + "RMG-database/"
conda_path = unc_folder + "conda/"
expt_yaml_file = ""
output_path  = "/scratch/blais.ch/methanol_results_2022_04_21/"

# make pdf pages obj
pp = PdfPages('report.pdf')

# get all objective function files
obj_func_files = glob.glob(os.path.join(output_path,"run_*/cantera/objective_function.txt"))
obj_func_files_log = glob.glob(os.path.join(output_path,"run_*/cantera/objective_function_log.txt"))

# run through all objective function files. also, return model size
obj_func_dict = {}
model_size_dict = {}
pattern = re.compile('run_+\d')
for file in obj_func_files:
    with open(file, "r") as f:
        line = f.readline()
        path,obj_func = line.split(":")
        
        # find the string for the run number 
        pattern = re.compile('run_\d{4}')
        match = re.search('run_\d{4}',file)
        run_num = int(match.group(0).replace('run_', ""))
        
        obj_func_dict[run_num] = float(obj_func)
        
    # access log file to get number fo species    
    new_addr = file.replace("cantera/objective_function.txt", "RMG.log")
    with open(new_addr, "r") as f:
        for line in f.readlines():
            if re.search('The final model core has', line, re.I):
                match_spec = re.search('[0-9]+ species',line)
                spec_num = int(match_spec.group(0).replace(" species", ""))
                
                match_reac = re.search('[0-9]+ reactions',line)
                reac_num = int(match_reac.group(0).replace(" reactions", ""))

                model_size_dict[run_num] = (spec_num,reac_num)
                
low_obj_func_dict = {}
low_model_size_dict = {}
for key, value in obj_func_dict.items(): 
    if value <= 1000.0 and value >= -1000.0:# or value >= 294:
        low_obj_func_dict[key] = value
        low_model_size_dict[key] = model_size_dict[key]

# run through all log objective function files
obj_func_dict_log = {}
model_size_dict_log = {}
for file in obj_func_files_log:
    with open(file, "r") as f:
        line = f.readline()
        path,obj_func = line.split(":")
        
        # find the string for the run number 
        pattern = re.compile('run_\d{4}')
        match = re.search('run_\d{4}',file)
        run_num = int(match.group(0).replace('run_', ""))
        
        obj_func_dict_log[run_num] = float(obj_func)
        
    # access log file to get number of species    
    new_addr = file.replace("cantera/objective_function_log.txt", "RMG.log")
    with open(new_addr, "r") as f:
        for line in f.readlines():
            if re.search('The final model core has', line, re.I):
                match_spec = re.search('[0-9]+ species',line)
                spec_num = int(match_spec.group(0).replace(" species", ""))
                
                match_reac = re.search('[0-9]+ reactions',line)
                reac_num = int(match_reac.group(0).replace(" reactions", ""))

                model_size_dict_log[run_num] = (spec_num,reac_num)
                
low_obj_func_dict_log = {}
low_model_size_dict_log = {}
for key, value in obj_func_dict_log.items(): 
    if value <= 1000.0 and value >= -1000.0:# or value >= 294:
        low_obj_func_dict_log[key] = value
        low_model_size_dict_log[key] = model_size_dict_log[key]


def plot_model_size(model_size_dict, obj_func_dict, pp, log=False, label=""):
    '''
    plot the species and reaction vs the objective function
    '''
    obj_list = []
    spec_list = []
    reac_list = []
    for expt, size in model_size_dict.items():
        obj_list.append(obj_func_dict[expt])
        spec_list.append(model_size_dict[expt][0])
        reac_list.append(model_size_dict[expt][1])

    plt.scatter(spec_list,obj_list)
    if log:
        plt.yscale('log')
    
    plt.title(f'{label} # of species vs objective function')
    pp.savefig(plt.gcf())
    plt.clf() 

    plt.scatter(reac_list,obj_list)
    if log:
        plt.yscale('log')

    plt.title(f'{label} # of reactions vs objective function')
    pp.savefig(plt.gcf())
    plt.clf() 

# plot em 
plot_model_size(
    model_size_dict, 
    obj_func_dict, 
    pp, 
    log=False,
    label="linear obj func")

plot_model_size(
    low_model_size_dict, 
    low_obj_func_dict, 
    pp, 
    log=False,
    label="linear obj func low")

plot_model_size(
    model_size_dict_log, 
    obj_func_dict_log, 
    pp, 
    log=False,
    label="log obj func")

plot_model_size(
    low_model_size_dict_log, 
    low_obj_func_dict_log, 
    pp, 
    log=False,
    label="log obj func low")

# load pickle file for perturbations
sobol_map = []
pickle_file = unc_folder + 'sobol_map.pickle' 
with open(pickle_file, "rb") as file:
    sobol_map = pickle.load(file)

sobol_map_copy = sobol_map 
for key, perts in sobol_map.items():
    sobol_map[key] = perts[1].tolist()
    
df = pd.DataFrame.from_dict(sobol_map)


# make this a function? for now will just plot them all as is

def plot_perturbs(df, obj_func_dict, pp, log=False, label=""):
    '''
    plot the perturbed values vs the objective function 
    '''
    counter = 0
    for col in df:
        counter += 1 
        obj_list = []
        val_list = []
        for run, obj in obj_func_dict.items():
            obj_list.append(obj)
            val_list.append(df[col][run])
            
        plt.scatter(val_list,obj_list)
        plt.title(f'{label} {col}')
        if log:
            plt.yscale('log') 
        pp.savefig(plt.gcf())
        plt.clf() 

# plot em
plot_perturbs(
    df, 
    obj_func_dict, 
    pp, 
    log=True,
    label="linear obj func")

plot_perturbs(
    df, 
    low_obj_func_dict, 
    pp, 
    log=True,
    label="linear obj func low")

plot_perturbs(
    df, 
    obj_func_dict_log, 
    pp, 
    log=False,
    label="log obj func")

plot_perturbs(
    df, 
    low_obj_func_dict_log, 
    pp, 
    log=False,
    label="log obj func low")
        


# generate list of passed and failed runs
# get number of folders
total_run_folders = glob.glob(os.path.join(output_path,"run_*/"))
all_runs = list(range(0,len(total_run_folders)))

passed_runs = list(obj_func_dict.keys())
passed_runs.sort()
passed_runs

# use a set to get the difference between passed and failed runs 
failed_runs_set = set(all_runs).difference(passed_runs)
failed_runs = list(failed_runs_set)
print("passed: ", len(passed_runs), "\nfailed: ", len(failed_runs))


# plot the passed vs failed runs 
counter = 0
for i,col in enumerate(df):
    counter += 1 
    failed_val_list = []
    passed_val_list = []
    for run in failed_runs:
        failed_val_list.append(df[col][run])
    for run in passed_runs:
        passed_val_list.append(df[col][run])

    plt.figure(i)   
    bw_adjust = 0.5
    sns.kdeplot(np.array(failed_val_list), bw_adjust=bw_adjust, label = "failed")
    sns.kdeplot(np.array(passed_val_list), bw_adjust=bw_adjust, label = "passed")
    plt.legend()
    plt.title(col)
    pp.savefig(plt.gcf())
    plt.clf() 


# get the values that are closest to the desired obj function value (0)
closest_linear = min(obj_func_dict, key=obj_func_dict.get)
print("closest linear obj func : {} : {closest_linear}")
# may not be the closest, we currently don't have any + values though
closest_log = max(obj_func_dict_log, key=obj_func_dict_log.get)


# # Analyze the model that fits best
run_num = str(closest_linear).zfill(4)
model_path =  os.path.join(output_path,f"run_{run_num}")
cti_path = os.path.join(model_path, "cantera/ct_analysis.csv")
data_graaf = pd.read_csv(cti_path)

ax = data_graaf.plot.scatter("graaf H2O TOF 1/s","RMG H2O TOF 1/s",loglog=True,color="b")
data_graaf.plot.scatter("graaf MeOH TOF 1/s","RMG MeOH TOF 1/s",loglog=True,color="r", ax=ax)
ax.set_ylim(1e-4, 5e-1)
ax.set_xlim(1e-4, 5e-1)

ax.set_xlabel("experimental TOF (1/s)")
ax.set_ylabel("simulated TOF (1/s)")
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.legend(labels=["MeOH TOF","H2O TOF"])

# make x=y line
plt.gca().set_xlim(left=1e-5)
plt.gca().set_ylim(bottom=1e-5)
xpoints = ypoints = plt.xlim()
plt.plot(xpoints, ypoints, linestyle='--', color='k', lw=3, scalex=False, scaley=False)

pp.savefig(plt.gcf())
plt.clf() 

# plot log obj function that's closest
run_num = str(closest_log).zfill(4)
model_path =  os.path.join(output_path,f"run_{run_num}")
cti_path = os.path.join(model_path, "cantera/ct_analysis.csv")

ax = data_graaf.plot.scatter("graaf H2O TOF 1/s","RMG H2O TOF 1/s",loglog=True,color="b")
data_graaf.plot.scatter("graaf MeOH TOF 1/s","RMG MeOH TOF 1/s",loglog=True,color="r", ax=ax)
ax.set_ylim(1e-4, 5e-1)
ax.set_xlim(1e-4, 5e-1)

ax.set_xlabel("experimental TOF (1/s)")
ax.set_ylabel("simulated TOF (1/s)")
ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.legend(labels=["MeOH TOF","H2O TOF"])


# make x=y line
plt.gca().set_xlim(left=1e-5)
plt.gca().set_ylim(bottom=1e-5)
xpoints = ypoints = plt.xlim()
plt.plot(xpoints, ypoints, linestyle='--', color='k', lw=3, scalex=False, scaley=False)
pp.savefig(plt.gcf())
plt.clf() 

pp.close()
















