import numpy as np
import yaml
import sys
import logging
import random
import time
sys.path.append(
    "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/")
from rmg_gua.gua_cantera.Spinning_basket_reactor.sbr import MinSBR

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


# def change_model(path, rdict):
#     """
#     change rms file to have new A and Ea values for reaction rxn_num
#     rdict = {"A_log_x or A_stick_x": A, "Ea_x": Ea} where x is rxn num
#     """

#     with open(path, "r") as f: 
#         ct_mech = yaml.load(f, Loader=yaml.FullLoader)

#     # get the most sensitive reaction from cti mech
#     for key, value in rdict.items():
#         num = int(key.split("_")[-1])
#         print(key, value, num)
#         if "A" in key and "log" in key:
#             A_i = 10**float(value)
#             ct_mech["surface1-reactions"][num]["rate-constant"]["A"] = A_i
#         elif "A" in key and "stick" in key:
#             A_i = float(value)
#             ct_mech["surface1-reactions"][num]["rate-constant"]["A"] = A_i
#         elif "E" in key:
#             Ea_i = float(value)
#             ct_mech["surface1-reactions"][num]["rate-constant"]["Ea"] = Ea_i
#         else:
#             logging.error(f"key {key} not recognized")

#     # save the yaml as a new file
#     new_path = path.replace(".yaml", "_new.yaml")
#     with open(new_path, "w") as f:
#         yaml.dump(ct_mech, f)

#     return new_path



#To use PEUQSE, you can have a function, but you also need to make a function wrapper that takes *only* the parameters as a single vector.
def simulationFunction(parameters):
    #here x is a scalar or an array and "a" and "b" are constants for the equation.
    """
    run rms reactor. 
    x = experimental data. I think we want it to run all of them at once? 
    logA = rxn 1 a-factor
    Ea = rxn 1 activation energy

    """
    t1 = time.asctime()
    t1s = time.time()
    # build and run the simulation
    file_path = "../../baseline/cantera/chem_annotated.cti"
    kin_par_path = "./ct_initial_small.yaml"
    expt_condts = "./expt_data.yaml"
    lookup_dict_file = "./rmg_2_ck_dict.yaml"
    CH3OH_X = []
    CO_X = []
    CO2_X = []
    H2_X = []
    H2O_X = []
    print("input parameters: ", parameters)
    # change the rms file. now doing all reactions in mechanism
    with open(kin_par_path, 'r') as file:
        kin_par_dat= yaml.safe_load(file)

    # unpack the peuquse parameters for use in modifying sim input
    input_a_ea = dict(
        zip(kin_par_dat["labels"], [parameters[0], parameters[1]]))


    with open(expt_condts, 'r') as file:
        data = yaml.safe_load(file)

    # pick just one experiment for example. can in the future use multiprocessing to solve faster, 
    # for now just do in series. 
    # get number of experiments: 
    run_test = True
    if run_test: 
        n_expts = 3
    else:
        n_expts = len(data["catalyst_area"])
    
    print(f"length is {n_expts} in sim")
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
              
            with open(lookup_dict_file, "r") as f: 
                lookup_dict =yaml.load(f, Loader = yaml.FullLoader)
            
        print(f"starting sim {run}")
        
        # we modify rates in memory
        test_sbr = MinSBR(
            file_path,
            reac_config=conditions,
            rtol=1.0e-11,
            atol=1.0e-22,
            new_rate_dict = input_a_ea
        )
        results = test_sbr.run_reactor_ss_memory()
        
        CH3OH_X.append(results[lookup_dict["CH3OH"]])
        CO_X.append(results[lookup_dict["CO"]])
        CO2_X.append(results[lookup_dict["CO2"]])
        H2_X.append(results[lookup_dict["H2"]])
        H2O_X.append(results[lookup_dict["H2O"]])
        
        


    y_data = np.vstack([CH3OH_X, CO_X, CO2_X, H2_X, H2O_X])
    
    print("Results:\n")
    print("CH3OH: ", CH3OH_X)
    print("CO: ", CO_X)
    print("CO2: ", CO2_X)
    print("H2: ", H2_X)
    print("H2O: ", H2O_X)
    print("sim done")
    t2 = time.asctime()
    t2s = time.time()
    
    print(f"start time: {t1}")
    print(f"end time: {t2}")
    print(f"elapsed: {t2s-t1s} seconds")
    return y_data 

def simulation_function_wrapper(parametersArray):#this has a and b in it.
    y = simulationFunction(parametersArray)  #an alternatie simpler syntax to unpack the parameters would be: simulationFunction(x_values_for_data, *parametersArray) 
    return y
