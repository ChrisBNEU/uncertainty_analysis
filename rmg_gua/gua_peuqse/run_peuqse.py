
import numpy as np 
import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput
import sys
sys.path.insert(0, '../../')



# Just a simple example. The user can also put the values in directly into the runfile or extract from a csv, for example.
# import observed_values_00
import rmg_model  # function for peuquse to optimize
import yaml


# load the experimental data yaml so we get a dict of lists for each parameter
expt_data = rmg_model.repackage_yaml("./rms_model/small_expt.yaml")
expt_data


# build the x-data array
x_data = []
x_data.append(expt_data["pressure"])
x_data.append(expt_data["species_CO"])
x_data.append(expt_data["species_CO2"])
x_data.append(expt_data["species_H2"])
x_data.append(expt_data["temperature"])
x_data.append(expt_data["volume_flowrate"])
x_data = np.array(x_data)



# build y-data array
y_data = []
y_data.append(expt_data["output_CH3OH"])
y_data.append(expt_data["output_H2O"])
y_data = np.array(y_data)
y_data


#Provide the observed X values and Y values and uncertainties -- all should be arrays or lists with nesting like [[1,2,3]] or [[1,2,3,4],[4,5,6,6]]
UserInput.responses['responses_abscissa'] = x_data
UserInput.responses['responses_observed'] = y_data

# don't know Graaf uncertainties currently, I think this is optional
# UserInput.responses['responses_observed_uncertainties'] = observed_values_00.observed_data_y_values_uncertainties



## initial values for As and Eas
import math

param_prior_values = {
    'A1': 8.36e+17, 'Ea1': 19479.16336492844,
    'A2': 4.17e+17, 'Ea2': 13071.514023053496,
    'A3': 4.18e+17, 'Ea3': 108372.25737624406,
    'A4': 4.18e+17, 'Ea4': 62236.621421176926,
    'A5': 4.18e+17, 'Ea5': 65143.36068000439,
}
param_prior_values_log = {
    'logA1': math.log10(param_prior_values['A1']), 'Ea1': param_prior_values['Ea1'],
    'logA2': math.log10(param_prior_values['A2']), 'Ea2': param_prior_values['Ea2'],
    'logA3': math.log10(param_prior_values['A3']), 'Ea3': param_prior_values['Ea3'],
    'logA4': math.log10(param_prior_values['A4']), 'Ea4': param_prior_values['Ea4'],
    'logA5': math.log10(param_prior_values['A5']), 'Ea5': param_prior_values['Ea5'],
}

param_log_unc = {
    'A1': 1.0, 'Ea1': 10000,
    'A2': 1.0, 'Ea2': 10000,
    'A3': 1.0, 'Ea3': 10000,
    'A4': 1.0, 'Ea4': 10000,
    'A5': 1.0, 'Ea5': 10000,
}

# convert to lists
param_prior_values_list = [val for val in param_prior_values.values()]
param_prior_values_log_list = [val for val in param_prior_values_log.values()]
param_log_unc_list = [val for val in param_log_unc.values()]

param_name_dict = {key:key for key in param_prior_values.keys()}


#Optional: provide labels for the responses axes and parameter names.
UserInput.simulated_response_plot_settings['x_label'] = ["pressure (Pa)", "mole frac CO", "mole frac CO2", "mole frac H2", "temperature (K)", "volume flowrate (m^3/s)"]
UserInput.simulated_response_plot_settings['y_label'] = ["CH3OH TOF (1/s)", "H2O TOF (1/s)"]
UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = param_name_dict

#Provide the prior distribution and uncertainties of the individual parameters.
# prior expected values for a and b
UserInput.model['InputParameterPriorValues'] = param_prior_values_log_list
# required. #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
UserInput.model['InputParametersPriorValuesUncertainties'] = param_log_unc_list
# UserInput.model['InputParameterInitialGuess'] = [150,400] #Can optionally change the initial guess to be different from prior means.
# UserInput.model['InputParameterPriorValues_upperBounds'] = [1, 20000] 
# UserInput.model['InputParameterPriorValues_lowerBounds'] = [0, 0]


#Provide a function that returns simulated values -- must of the same form as observed values, should be arrays or lists with nesting like [[1,2,3]] or [[1,2,3,4],[4,5,6,6]]
# This must simulate with *only* the parameters listed above, and no other arguments.
UserInput.model['simulateByInputParametersOnlyFunction'] = rmg_model.simulation_function_wrapper

#mcmc length should typically be on the order of 10,000 per parameter. By default, the burn in will be the first 10% of the mcmc length.
# 10000 is the default.
UserInput.parameter_estimation_settings['mcmc_length'] = len(param_prior_values_log_list)*10000

#After filinlg the variables of the UserInput, now we make a 'parameter_estimation' object from it.
PE_object = PEUQSE.parameter_estimation(UserInput)

#Now we can do the mcmc!
PE_object.doMetropolisHastings()
#Another option would be PE_object.doEnsembleSliceSampling(), one can also do grid search or an astroidal distribution search.

#Finally, create all plots!
PE_object.createAllPlots()
#The createAllPlots function calls each of the below functions so that the user does not have to.
#    PE_object.makeHistogramsForEachParameter()
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlots()

#########Optional example of saving and loading PE_objects after running the mcmc.
#########This feature requires having dill installed (pip install dill, https://pypi.org/project/dill/)
try:
    import dill
    dillModuleExists = True
except:
    dillModuleExists = False

#Optionally, one can save a PE_object for later, if the dill module has been installed.
if dillModuleExists == True:
    PE_object.save_to_dill("PE_object_00a0")
    #to load a PE_object after some time, first one has to put (any) UserInput to create a PE_object, then to load from file.

    #Normally, we would do the loading and plotting in another python file, but for this example the syntax is being demonstrated below within the same file.
    PE_object2 = PEUQSE.parameter_estimation(UserInput)
    PE_object2 = PE_object2.load_from_dill("PE_object_00a0")
    PE_object2.createAllPlots()

 
# # Misc functions for creating inputs for PEUQSE


listy1 = [1,2,3,4,5]
listy2 = [7,8,9,10,11]

np.vstack([listy1, listy2])


import os
import shutil
import sys
import yaml



def change_model(path, A, Ea, rxn_num=0):
    """
    change rms file to have new A and Ea values for reaction rxn_num
    """
    
    with open(path, "r") as f:
         rms_mech = yaml.load(f, Loader=yaml.FullLoader)

    rms_mech["Reactions"][rxn_num]["kinetics"]["A"] = A
    rms_mech["Reactions"][rxn_num]["kinetics"]["Ea"] = Ea

    # save the yaml as a new file
    new_path = path.replace(".rms", "_new.rms")
    with open(new_path, "w") as f:
        yaml.dump(rms_mech, f)




A = 0.5
Ea = 1
rxn_num = 1
path = "./rms_model/rms/chem25.rms"
new_path = path.replace(".rms", "_new.rms")
rms_mech = yaml.load(open(path, "r"), Loader=yaml.FullLoader)
# print ("old mech: \n", rms_mech["Reactions"][rxn_num]["kinetics"], "\n")
# change_model(path, A, Ea, rxn_num=rxn_num)

# rms_mech = yaml.load(open(new_path, "r"), Loader=yaml.FullLoader)
# print("new mech: \n", rms_mech["Reactions"][rxn_num]["kinetics"], "\n")







# take the most sensitive reactions and perturb their A's 
import pickle
sens_pickle = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/" + \
    "uncertainty_analysis/rmg_gua/baseline/sens_rms_dict.pickle"

# do the top 5 sensitive rxns
with open(sens_pickle, "rb") as f:
    # load pickle file
    sens_dict = pickle.load(f)
counter = 0

most_sens_list = {}
for key, value in sens_dict.items(): 
    counter +=1
    if counter > 5:
        break
    print(key, value[0])
    most_sens_reac = key.split("<=>", )
    reac = key.split("<=>")[0].split("+")
    prod = key.split("<=>")[1].split("+")
    most_sens_list[key] = [reac, prod]
most_sens_list


for key, val in most_sens_list.items(): 
    print(len(val))
    val


type(rms_mech)


# get the rxn_num for each reaction in the 
import copy
most_sens_list_copy = copy.deepcopy(most_sens_list)
for num, rxn in enumerate(rms_mech["Reactions"]):
    for reaction, specs in most_sens_list.items():
        if set(specs[0]) == set(rxn["reactants"]) and set(specs[1]) == set(rxn["products"]):
            if len(most_sens_list_copy[reaction]) >=3:
                print("more than 1 match for ", reaction, num)
                pass
            else:
                # print(num, rxn["reactants"], rxn["products"], rxn["kinetics"])
                most_sens_list_copy[reaction].append(num)

pert_rxns = [val[2] for val in most_sens_list_copy.values()]
pert_rxns



for num in pert_rxns:
    print(rms_mech["Reactions"][num]["kinetics"], rms_mech["Reactions"][num]["reactants"], rms_mech["Reactions"][num]["products"])


expt_condts = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/gua_cantera/all_experiments_reorg_sbr.yaml"
with open(expt_condts, "r") as f:
    expt_list = yaml.load(f, Loader=yaml.FullLoader)


import math
rtol = 0.1
def compare_dicts(dict1, dict2):
    match = (
        dict1["use_for_opt"] == dict2["use_for_opt"],

        math.isclose(dict1["catalyst_area"], dict2["catalyst_area"], rel_tol=rtol, abs_tol=1e-08,),
        math.isclose(dict1["pressure"],dict2["pressure"],rel_tol=rtol, abs_tol=1e-08,),
        math.isclose(dict1["temperature"],dict2["temperature"],rel_tol=rtol, abs_tol=1e-08,),
        math.isclose(dict1["volume"], dict2["volume"], rel_tol=rtol, abs_tol=1e-08,),
        math.isclose(dict1["volume_flowrate"],dict2["volume_flowrate"],rel_tol=rtol, abs_tol=1e-08,),
        math.isclose(
            dict1["species"]["CO"],
            dict2["species"]["CO"],rel_tol=rtol, abs_tol=1e-08,),
        math.isclose(
            dict1["species"]["CO2"],
            dict2["species"]["CO2"], rel_tol=rtol, abs_tol=1e-08,),
        math.isclose(
            dict1["species"]["H2"],
            dict2["species"]["H2"], rel_tol = rtol, abs_tol = 1e-08,),
        )
    if all(match):
        return True
        print("match")
    else:
        return False


matches = {}
for num,exp in enumerate(expt_list[:212]): 
    for next_num, next_exp in enumerate(expt_list[num+1:]): 
        # if next_num < 5:
        #     print("next_expt:" , next_num)
        match_exp = compare_dicts(exp, next_exp)
        if match_exp:
            if num in matches.keys():
                matches[num].append(num + next_num)
            else:
                matches[num] = [num + next_num]
max_len = 0
match_expt = None
for expt, matchy in matches.items(): 
    if len(matchy) > max_len:
        max_len = len(matchy)
        match_expt = expt

print("max_len", max_len)
print("match_expt", match_expt)
print("matches", matches[match_expt])
    

new_expt_yaml = []

new_expt_yaml.append(expt_list[match_expt])
for expt in matches[match_expt]:
    new_expt_yaml.append(expt_list[expt])

new_expt_file = "./rms_model/small_expt.yaml"
with open(new_expt_file, "w") as f:
    yaml.dump(new_expt_yaml, f)


match_exp = compare_dicts(expt_list[0], expt_list[0])
match_exp


matches


expt_list[0]





