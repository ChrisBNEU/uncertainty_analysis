#!/usr/bin/env python
# coding: utf-8

# ## get all of the runs for global uncertainty
# 
# UQTK requires 3 inputs: 
# - a text file containing all of the perturbed values
# - a text file containing all of the ranges for the perturbed values
# - a text file containing the outputs that you'd like to get sensitivity to (for us, meoh tof, h2o tof)
# 
# There is a split for training vs validation, I am just goint to do 75% to start idk what best practices are

import pandas as pd
import csv
import pickle
import os
import copy

import numpy as np
import math
import csv
import glob
import re

def create_files(
    output_path,
    M,
    ):
    
    # # metal BE values from RMG-DB
    # output_path  = "/scratch/blais.ch/methanol_results_2022_08_21_reduced/"

    # we need to know which runs failed. use script from the postprocessing notebook that 
    # checks for the objective function file
    # we only want the ct_analysis files that were completed most recently, 
    # which will have an objective function file with them
    obj_func_filename = "objective_function_log_ct_analysis.txt"
    ct_filename = "ct_analysis.csv"
    ct_data_files = glob.glob(os.path.join(output_path,f"run_*/cantera/{ct_filename}"))

    obj_func_files = [i.replace(ct_filename, obj_func_filename) for i in ct_data_files]
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
        new_addr = file.replace(f"cantera/{obj_func_filename}", "RMG.log")
        with open(new_addr, "r") as f:
            for line in f.readlines():
                if re.search('The final model core has', line, re.I):
                    match_spec = re.search('[0-9]+ species',line)
                    spec_num = int(match_spec.group(0).replace(" species", ""))

                    match_reac = re.search('[0-9]+ reactions',line)
                    reac_num = int(match_reac.group(0).replace(" reactions", ""))

                    model_size_dict[run_num] = (spec_num,reac_num)

    # generate list of passed and failed runs
    # get number of folders

    total_run_folders = glob.glob(os.path.join(output_path,"run_*/"))
    all_runs = list(range(0,len(total_run_folders)))

    passed_runs = list(obj_func_dict.keys())
    passed_runs.sort()

    print("number of passed runs: ", len(passed_runs))

    # use a set to get the difference between passed and failed runs 
    failed_runs_set = set(all_runs).difference(passed_runs)
    failed_runs = list(failed_runs_set)

    closest = min(obj_func_dict, key=lambda y: abs(obj_func_dict[y]))
    print("run with closest match: ", closest, obj_func_dict[closest])

    # # sanitize inputs
    # first, load in the perturbed values  
    # - the sobol map pickle has the value of the perturbed parameters. load this, then scale the a factors so we're looking at their exponent? 
    # - everything else should be okay I think? 

    sobol_path = f"{output_path}sobol_map.pickle"
    with open(sobol_path, "rb") as input_file:
        perturb_values = pickle.load(input_file)

    # reconstruct dict so it is {key:list of values}
    pert_val_san = {}
    for param, values in perturb_values.items():
        # extract a list of just perturbed values, since we 
        # have a tuple with that and other metadata
        just_pert = values[1].tolist()

        # check for A factors that are not sticking coefficient
        # record the value at
        if min(just_pert) > 1.0 and param.endswith("/A"): 
            just_pert = [math.log(i,10) for i in just_pert]

        pert_val_san[param] = just_pert

    # re-organize so we have a list of dicts, with the perturbed values as each individual dict? 
    pert_val_san_list = []
    for i in range(M):
        ind_pert_dict = {}
        for key in pert_val_san.keys():
            ind_pert_dict[key] = pert_val_san[key][i]

        pert_val_san_list.append(ind_pert_dict)

    ## passed runs
    pert_val_san_list_passed = []
    for run, data in enumerate(pert_val_san_list):
        if run in passed_runs:
            pert_val_san_list_passed.append(data)

    # load the range pickle file
    sobol_range_path = f"{output_path}sobol_range_map.pickle"
    with open(sobol_range_path, "rb") as input_file:
        perturb_range_values = pickle.load(input_file)

    # now just make a range text file. 
    range_file = "parameter_ranges_meoh.txt"
    with open(range_file, "w") as f:
        for key, value in perturb_range_values.items():
            f.write("{0:f} {1:f}\n".format(value[1], value[2]))


    # # get parameter names
    pert_value_names = list(perturb_range_values.keys())

    # write to output file 
    with open("parameter_names_meoh.txt", "w") as f:
        for name in pert_value_names:
            f.write(name+"\n")


    # # get output data and write to file
    ct_data_dict = {}

    data_keys = [
        "RMG MeOH TOF 1/s",
        "RMG H2O TOF 1/s",
    ]

    # we have a number of completed runs, but we really need the number 
    # of runs where meoh/h2o rop are greaterthan 0. 

    for file in ct_data_files:
        # find the string for the run number 
        pattern = re.compile('run_\d{4}')
        match = re.search('run_\d{4}',file)
        run_num = int(match.group(0).replace('run_', ""))
        df = pd.read_csv(file)
        data_run = df[(df['Unnamed: 0'] == 80)]
        data_dict = {}

        for key in data_keys: 
            output_value = float(data_run[key])
            if  output_value <=0: 
                output_value = 0
            else: 
                output_value = math.log(output_value, 10)
            data_dict[key] = output_value


        ct_data_dict[run_num]= data_dict

    # sort by run

    # make dict of passed runs with their perturbed values
    pert_val_san_dict_passed = {}
    for run, data in enumerate(pert_val_san_list):
        if run in passed_runs:
            pert_val_san_dict_passed[run] = data

    # match up all of the output data for the cantera runs that passed
    pert_val_san_dict_passed_copy = pert_val_san_dict_passed.copy()
    ct_data_dict_copy = ct_data_dict.copy()
    for run, data in ct_data_dict.items():

        for key, value in ct_data_dict_copy[run].items():
            if value == 0:
                del pert_val_san_dict_passed_copy[run]
                del ct_data_dict_copy[run]
                print("deleted run: ", run, " because it has TOF of 0")
                break

    # now just make a text file. 
    output_file = "outputs_meoh.txt"
    with open(output_file, "w") as f:
        for key, data_keys in ct_data_dict_copy.items():
            val_str_list = []
            for value in data_keys:
                val_str_list.append("{0:e}".format(data_keys[value]))
            value_str = " ".join(val_str_list)
            f.write(f"{value_str}\n")

    # make the input file now that we know which runs to keep
    csv_file = "Input_meoh.txt"
    csv_columns = list(perturb_values.keys())



    try:
        with open(csv_file, 'w') as f:
            writer = csv.writer(f, delimiter=" ")
            # we will remove the header row, seems dicy because I do not know if the 
            # values are written correctly
            # writer.writeheader()
            for key,perts in pert_val_san_dict_passed_copy.items():
                data = []
                for value in perts.values():
                    data.append(value)
                writer.writerow(data)

    except IOError:
        print("I/O error")


    # update workflow.x file with training and valudation numbers
    num_passed = len(pert_val_san_dict_passed)
    num_train = round(0.75*len(passed_runs))
    num_val = num_passed-num_train
    print("Training: ", num_train,"\nValidation: ", num_val)


    with open('workflow_meoh.x', 'r') as f:
        # read a list of lines into data
        data = f.readlines()

    for line, datum in enumerate(data): 
        if "TRAIN=" in datum:
            data[line] = f"TRAIN={num_train}\n"
        if ("VAL=" in datum) and not ("NVAL=" in datum):
            data[line] = f"VAL={num_val}\n"      

    with open('workflow_meoh.x', 'w') as f:
        f.writelines(data)
        

if __name__ == "__main__":
    output_path  = "/scratch/blais.ch/methanol_results_2022_08_25_just_thermo_5000/"
    M = 5000
    create_files(
    output_path,
    M=M,
    )