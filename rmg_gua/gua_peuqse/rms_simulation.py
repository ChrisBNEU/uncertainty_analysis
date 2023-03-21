import numpy as np
import yaml
import sys
import logging
import random
sys.path.append(
    "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/")
from rmg_gua.gua_rms.sbr import rms_sbr

def repackage_yaml(yaml_file):
    """
    make the yaml for expt data dict with lists for entries
    """
    with open(yaml_file, "r") as f:
        expt_yaml = yaml.load(f, Loader=yaml.FullLoader)
    
    # get the keys from the first entry. 
    new_expt_dict = {}
    for num, expt in enumerate(expt_yaml): 
        for key, value in expt.items():
            if not isinstance(value, dict): 
                if key not in new_expt_dict.keys():
                    new_expt_dict[key] = [value]
                else: 
                    new_expt_dict[key].append(value)
            else: 
                # bump out nested dicts to the top level (e.g. concentrations in)
                for sub_key, sub_value in value.items(): 
                    new_key = key + "_" + sub_key
                    if new_key not in new_expt_dict.keys():
                        new_expt_dict[new_key] = [sub_value]
                    else:
                        new_expt_dict[new_key].append(sub_value)
                    
    return new_expt_dict
    # make the 


def change_model(path, rdict):
    """
    change rms file to have new A and Ea values for reaction rxn_num
    rdict = {"A_log_x or A_stick_x": A, "Ea_x": Ea} where x is rxn num
    """

    with open(path, "r") as f:
        rms_mech = yaml.load(f, Loader=yaml.FullLoader)

    for key, value in rdict.items():
        num = int(key.split("_")[-1])

        if "A" in key and "log" in key:
            A_i = 10**float(value)
            rms_mech["Reactions"][num]["kinetics"]["A"] = A_i
        elif "A" in key and "stick" in key:
            A_i = float(value)
            rms_mech["Reactions"][num]["kinetics"]["A"] = A_i
        elif "E" in key:
            Ea_i = float(value)
            rms_mech["Reactions"][num]["kinetics"]["Ea"] = Ea_i
        else:
            logging.error(f"key {key} not recognized")

    # save the yaml as a new file
    new_path = path.replace(".rms", "_new.rms")
    with open(new_path, "w") as f:
        yaml.dump(rms_mech, f)

    return new_path



#To use PEUQSE, you can have a function, but you also need to make a function wrapper that takes *only* the parameters as a single vector.
def simulationFunction(parameters):
    #here x is a scalar or an array and "a" and "b" are constants for the equation.
    """
    run rms reactor. 
    x = experimental data. I think we want it to run all of them at once? 
    logA = rxn 1 a-factor
    Ea = rxn 1 activation energy

    """
    # build and run the simulation
    file_path = "../baseline/rms/chem53.rms"
    kin_par_path = "./rms_initial_small.yaml"
    rmg_db_folder = "/Users/blais.ch/Documents/_01_code/RMG_env_1/RMG-database/"
    expt_condts = "expt_data.yaml"
    CH3OH_X = []
    CO_X = []
    CO2_X = []
    H2_X = []
    H2O_X = []
    # change the rms file. now doing all reactions in mechanism
    with open(kin_par_path, 'r') as file:
        kin_par_dat= yaml.safe_load(file)

    with open(file_path, 'r') as file: 
        rms_mech = yaml.safe_load(file)

    n_reactions = len(rms_mech["Reactions"])

    # unpack the peuquse parameters for use in modifying sim input
    input_a_ea = dict(
        zip(kin_par_dat["labels"], kin_par_dat["initial_values"]))

    file_path_new = change_model(file_path, input_a_ea)

    with open(expt_condts, 'r') as file:
        data = yaml.safe_load(file)

    # pick just one experiment for example. can in the future use multiprocessing to solve faster, 
    # for now just do in series. 
    # get number of experiments: 
    n_expts = len(data["catalyst_area"])
    for run in range(0,n_expts): 
        conditions = {}
        for key in data.keys():
            if "species_" in key and not "out" in key:
                spec_name = key.split("_")[-1]
                if "species" not in conditions.keys():
                    spec_name = key.split("_")[-1]
                    conditions["species"] = {}
                    conditions["species"][spec_name] = data[key][run]
                else:
                    conditions["species"][spec_name] = data[key][run]
            elif "output_" in key:
                spec_name = key.split("_")[-1]
                if "output" not in conditions.keys():
                    spec_name = key.split("_")[-1]
                    conditions["output"] = {}
                    conditions["output"][spec_name] = data[key][run]
                else:
                    conditions["output"][spec_name] = data[key][run]
                
            elif "species_out_" in key:
                spec_name = key.split("_")[-1]
                if "species_out" not in conditions.keys():
                    spec_name = key.split("_")[-1]
                    conditions["species_out"] = {}
                    conditions["species_out"][spec_name] = data[key][run]
                else:
                    conditions["species_out"][spec_name] = data[key][run]
            else: 
                conditions[key] = data[key][run]
              
        # print(conditions)
        print(f"starting sim {run}")
        test_sbr = rms_sbr(
            file_path_new,
            reac_config=conditions,
            rtol=1.0e-11,
            atol=1.0e-22,
        )
        results = test_sbr.run_simulation()
        
        CH3OH_X.append(results["CH3OH"])
        CO_X.append(results["CO"])
        CO2_X.append(results["CO2"])
        H2_X.append(results["H2"])
        H2O_X.append(results["H2O"])
        results = test_sbr.run_simulation()
        
        # CH3OH_X.append(random.random())
        # CO_X.append(random.random())
        # CO2_X.append(random.random())
        # H2_X.append(random.random())
        # H2O_X.append(random.random())

    y_data = np.vstack([CH3OH_X, CO_X, CO2_X, H2_X, H2O_X])
    print("sim done")
    return y_data 

def simulation_function_wrapper(parametersArray):#this has a and b in it.
    y = simulationFunction(parametersArray)  #an alternatie simpler syntax to unpack the parameters would be: simulationFunction(x_values_for_data, *parametersArray) 
    return y