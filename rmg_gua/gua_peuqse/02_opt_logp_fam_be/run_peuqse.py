import copy
import pickle
import shutil
import os
import math
import yaml

import numpy as np
import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput
import sys
repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(repo_dir)
results_path = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "02_opt_logp_fam_be", "config")
from rmg_gua.gua_peuqse.runtime_utilities import get_all_param_lists
run_test = False


import ct_simulation  # function for peuquse to optimize

if __name__ == "__main__":
    # Just a simple example. The user can also put the values in directly into the runfile or extract from a csv, for example.
    print("running job")

    results_path = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "02_opt_logp_fam_be", "config")
    # get the peuqse parameters
    param_dict = get_all_param_lists(results_path=results_path)     

    # load the exp data yaml
    with open("expt_data.yaml", "r") as f:
        expt_data = yaml.load(f, Loader=yaml.FullLoader)
        
    # load uncertainty yaml
    with open("expt_unc.yaml", "r") as f:
        expt_unc = yaml.load(f, Loader=yaml.FullLoader)
    
    x_data = []
    for expt in range(len(expt_data["catalyst_area"])):
        row = []
        row.append(expt_data["catalyst_area"][expt])
        row.append(expt_data["pressure"][expt])
        row.append(expt_data["species_CO"][expt])
        row.append(expt_data["species_CO2"][expt])
        row.append(expt_data["species_H2"][expt])
        row.append(expt_data["temperature"][expt])
        row.append(expt_data["volume_flowrate"][expt])
        x_data.append(row)
    x_data = np.array(x_data)
    print(f"length is {len(x_data[0])} in main")

    # build y-data array
    y_data = []
    y_data.append(expt_data['species_out_CH3OH'])
    y_data.append(expt_data['species_out_CO'])
    y_data.append(expt_data['species_out_CO2'])
    y_data.append(expt_data['species_out_H2'])
    y_data.append(expt_data['species_out_H2O'])
    y_data = np.array(y_data)


    # build y uncertainties array
    y_unc = []
    y_unc.append(expt_unc['species_out_CH3OH'])
    y_unc.append(expt_unc['species_out_CO'])
    y_unc.append(expt_unc['species_out_CO2'])
    y_unc.append(expt_unc['species_out_H2'])
    y_unc.append(expt_unc['species_out_H2O'])
    y_unc = np.array(y_unc)

    # Provide the observed X values and Y values and uncertainties -- all should be arrays or lists with nesting like [[1,2,3]] or [[1,2,3,4],[4,5,6,6]]
    UserInput.responses['responses_abscissa'] = x_data
    UserInput.responses['responses_observed'] = y_data

    # estimated graaf uncertainties. 
    UserInput.responses['responses_observed_uncertainties'] = y_unc

    #Optional: provide labels for the responses axes and parameter names.
    UserInput.simulated_response_plot_settings['x_label'] = [
        "catalyst_area", "pressure (Pa)", "mole frac CO", "mole frac CO2", "mole frac H2", "temperature (K)", "volume flowrate (m^3/s)"]
    UserInput.simulated_response_plot_settings['y_label'] = [
        "CH3OH mole Frac out", "CO mole Frac out", "CO2 mole Frac out", "H2 mole Frac out", "H2O mole Frac out"]
    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = param_dict["label_list"]

    #Provide the prior distribution and uncertainties of the individual parameters.
    # prior expected values for a and b
    UserInput.model['InputParameterPriorValues'] = param_dict["value_list"]
    # required. #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    UserInput.model['InputParametersPriorValuesUncertainties'] = param_dict["unc_list"]
    # UserInput.model['InputParameterInitialGuess'] = #Can optionally change the initial guess to be different from prior means.
    UserInput.model['InputParameterPriorValues_lowerBounds'] = param_dict["lower_list"]
    UserInput.model['InputParameterPriorValues_upperBounds'] = param_dict["upper_list"]
    
    UserInput.model['simulateByInputParametersOnlyFunction'] = ct_simulation.simulation_function_wrapper

    UserInput.parameter_estimation_settings['mcmc_threshold_filter_samples'] = True

    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0
    UserInput.parameter_estimation_settings['multistart_searchType'] = 'doOptimizeLogP'
    UserInput.parameter_estimation_settings['multistart_passThroughArgs'] = {'method':'BFGS'} #Here, BFGS is used. However, Nelder-Mead is usually what is recommended.
    UserInput.parameter_estimation_settings['multistart_initialPointsDistributionType'] = 'grid'
    UserInput.parameter_estimation_settings['multistart_exportLog'] = True
    
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)
    
    # run the model
    PE_object.doMultiStart()
    PE_object.createAllPlots() #This function calls each of the below functions so that the user does not have to.


