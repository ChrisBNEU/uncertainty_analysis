import numpy as np
import yaml
import sys
sys.path.append(
    "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/")
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


def change_model(path, logA, Ea, rxn_num=0):
    """
    change rms file to have new A and Ea values for reaction rxn_num
    """

    with open(path, "r") as f:
        rms_mech = yaml.load(f, Loader=yaml.FullLoader)

    if all([isinstance(logA, list), isinstance(Ea, list), isinstance(rxn_num, list)]):
        for num, rxn_num in enumerate(rxn_num):
            A_i = 10**float(logA[num])
            Ea_i = float(Ea[num])
            rms_mech["Reactions"][rxn_num]["kinetics"]["A"] = A_i
            print(f"A{num}", A_i)
            rms_mech["Reactions"][rxn_num]["kinetics"]["Ea"] = Ea_i
            print(f"Ea{num}", Ea_i)
    else: 
        A = 10**float(logA)
        Ea = float(Ea)
        rms_mech["Reactions"][rxn_num]["kinetics"]["A"] = A
        rms_mech["Reactions"][rxn_num]["kinetics"]["Ea"] = Ea

    # save the yaml as a new file
    new_path = path.replace(".rms", "_new.rms")
    with open(new_path, "w") as f:
        yaml.dump(rms_mech, f)

    return new_path



#To use PEUQSE, you can have a function, but you also need to make a function wrapper that takes *only* the parameters as a single vector.
def simulationFunction(logA1,Ea1,logA2,Ea2,logA3,Ea3,logA4,Ea4,logA5,Ea5):
    #here x is a scalar or an array and "a" and "b" are constants for the equation.
    """
    run rms reactor. 
    x = experimental data. I think we want it to run all of them at once? 
    logA = rxn 1 a-factor
    Ea = rxn 1 activation energy

    """
    # build and run the simulation
    file_path = "./rms_model/rms/chem25.rms"
    rmg_db_folder = "/Users/blais.ch/Documents/_01_code/RMG_env_1/RMG-database/"
    expt_condts = "./rms_model/small_expt.yaml"
    meoh_tof = []
    h2o_tof = []
    # change the rms file. rxn 21 is the CH3OH desporption rxn

    # manualy putting in most sensitive rxns in mechanism. 
    rxn_nums = [100, 106, 84, 38, 30]
    logA_list = [logA1, logA2, logA3, logA4, logA5]
    Ea_list = [Ea1, Ea2, Ea3, Ea4, Ea5]
    file_path_new = change_model(file_path, logA_list, Ea_list, rxn_num=rxn_nums)

    with open(expt_condts, 'r') as file:
        data = yaml.safe_load(file)

    # pick just one experiment for example. can in the future use multiprocessing to solve faster, 
    # for now just do in series. 
    for run, conditions in enumerate(data): 

        test_sbr = rms_sbr(
            file_path_new,
            reac_config=conditions,
            rtol=1.0e-11,
            atol=1.0e-22,
        )
        results = test_sbr.run_simulation()

        meoh_tof.append(results['RMG MeOH TOF 1/s'])
        h2o_tof.append(results['RMG H2O TOF 1/s'])
    
    y_data = np.vstack([meoh_tof, h2o_tof])

    return y_data

def simulation_function_wrapper(parametersArray):#this has a and b in it.
    print("parametersArray: ", parametersArray, "\nType : ", type(parametersArray))
    y = simulationFunction(*parametersArray)  #an alternatie simpler syntax to unpack the parameters would be: simulationFunction(x_values_for_data, *parametersArray) 
    print(y)
    return y