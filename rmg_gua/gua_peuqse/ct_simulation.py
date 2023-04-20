import numpy as np
import yaml
import sys
import os
import logging
import random
import time
import h5py
repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(repo_dir)

from rmg_gua.gua_cantera.Spinning_basket_reactor.sbr import MinSBR
from cantera import CanteraError
from PEUQSE import parallel_processing

def sim_init(project_path):
    global results_path
    global peuq_path
    global data
    global lookup_dict

    results_path = os.path.join(project_path, "config")
    peuq_path = os.path.join(project_path, "peuqse_out")

    expt_condts = os.path.join(results_path, "ct_expt_list.yaml")
    lookup_dict_file = os.path.join(results_path, "rmg_2_ck_dict.yaml")

    # load the exp data and lookup dict
    with open(expt_condts, 'r') as f:
        data = yaml.load(f, Loader = yaml.FullLoader)
    with open(lookup_dict_file, 'r') as f:
        lookup_dict = yaml.load(f, Loader = yaml.FullLoader)


# To use PEUQSE, you can have a function, but you also need to make a function wrapper that takes *only* the parameters as a single vector.
def simulationFunction(parameters, debug=False):
    #here x is a scalar or an array and "a" and "b" are constants for the equation.
    """
    run rms reactor. 
    """
    t1 = time.asctime()
    t1s = time.time()
    # build and run the simulation
    file_path = os.path.join(repo_dir, "rmg_gua", "baseline", "cantera", "chem_annotated.cti")
    

    CH3OH_X = []
    CO_X = []
    CO2_X = []
    H2_X = []
    H2O_X = []
    # change the rms file. now doing all reactions in mechanism


    # pick just one experiment for example. can in the future use multiprocessing to solve faster, 
    # for now just do in series. 
    # get number of experiments: 
    run_test = False
    n_expts = len(data)
            
    print(f"length is {n_expts} in sim")
    for run in range(0,n_expts): 
        conditions = data[run]
        print(f"starting sim {run}")
        
        # we modify rates in memory
        test_sbr = MinSBR(
            file_path,
            reac_config=conditions,
            rtol=1.0e-11,
            atol=1.0e-22,
            reaction_list=parameters, 
            results_path=results_path,
        )
        # try except loop for this because cantera errors out for certain A and Ea values
        try:
            results = test_sbr.run_reactor_ss_memory()
            CH3OH_X.append(results[lookup_dict["CH3OH"]])
            CO_X.append(results[lookup_dict["CO"]])
            CO2_X.append(results[lookup_dict["CO2"]])
            H2_X.append(results[lookup_dict["H2"]])
            H2O_X.append(results[lookup_dict["H2O"]])
        except Exception as e: #CanteraError: making bare exception just to see if there is something catchable 
            with open("./mpiproblem.txt", "w") as f:
                f.write(f"error on cantera run {run}, encountered exception {e}")
            print("could not solve system of equations with parameters [parameters], so setting outlet moles to nan")
            CH3OH_X.append(float('nan'))
            CO_X.append(float('nan'))
            CO2_X.append(float('nan'))
            H2_X.append(float('nan'))
            H2O_X.append(float('nan'))

    y_data = np.vstack([CH3OH_X, CO_X, CO2_X, H2_X, H2O_X])
    pnum = parallel_processing.currentProcessorNumber
    
    # is there a faster filetype? json?
    peuq_output_yaml_path = os.path.join(peuq_path, f"{pnum}_out.yaml")
    if os.path.exists(peuq_output_yaml_path):
        with open(peuq_output_yaml_path, "r") as f:
            peuq_output_yaml = yaml.safe_load(f)
    else: 
        peuq_output_yaml = []
        
    # not sure why but peuqse has an array for first iteration and a tuple for
    # successive ones for input parameters
    # print("dtype output: ", type(parameters[0]), type(y_data[0][0]))
    if isinstance(parameters, tuple):
        parameters = np.array(parameters).astype(float).tolist()
        y_data_list = y_data.astype(float).tolist()
    else: 
        parameters = parameters.astype(float).tolist()
        y_data_list = y_data.astype(float).tolist()
        
    
    peuq_output_yaml.append((parameters, y_data_list))
    
    with open(peuq_output_yaml_path, "w") as f:
        yaml.safe_dump(peuq_output_yaml, f)
    
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
    
    if debug:
        return test_sbr, y_data
    else:
        return y_data 
    

def simulation_function_wrapper(parametersArray):#this has a and b in it.
    y = simulationFunction(parametersArray)  #an alternate simpler syntax to unpack the parameters would be: simulationFunction(x_values_for_data, *parametersArray) 
    return y
