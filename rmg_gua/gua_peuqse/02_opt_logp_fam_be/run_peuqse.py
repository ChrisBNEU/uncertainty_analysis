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
from rmg_gua.gua_peuqse.runtime_utilities import get_all_param_lists, make_exp_data_lists
run_test = False


import ct_simulation  # function for peuquse to optimize

if __name__ == "__main__":
    # Just a simple example. The user can also put the values in directly into the runfile or extract from a csv, for example.
    print("running job")

    results_path = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "02_opt_logp_fam_be", "config")
    # get the peuqse parameters
    param_dict = get_all_param_lists(results_path=results_path)     

    # load the exp data yaml
    expt_data_path = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "02_opt_logp_fam_be", "expt_data.yaml")
    with open(expt_data_path, "r") as f:
        expt_data = yaml.load(f, Loader=yaml.FullLoader)
        
    # load uncertainty yaml. experimenting with different layout for X. 
    expt_unc_path = os.path.join(repo_dir, "rmg_gua", "gua_peuqse", "02_opt_logp_fam_be", "expt_unc.yaml")
    with open(expt_unc_path, "r") as f:
        expt_unc = yaml.load(f, Loader=yaml.FullLoader)
    
    data_path = os.path.join(repo_dir, "rmg_gua", "gua_cantera")
    _, x_data, y_data, y_unc = make_exp_data_lists(data_path, results_path, use_count=3)
    
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    y_unc = np.array(y_data)
    print(f"length is {len(x_data)} in main")

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
    # UserInput.model['InputParameterPriorValues_lowerBounds'] = param_dict["lower_list"]
    # UserInput.model['InputParameterPriorValues_upperBounds'] = param_dict["upper_list"]
    
    UserInput.model['simulateByInputParametersOnlyFunction'] = ct_simulation.simulation_function_wrapper

    UserInput.parameter_estimation_settings['mcmc_length'] = 10000
    
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0
    
    UserInput.parameter_estimation_settings['multistart_searchType'] = 'doEnsembleSliceSampling'
    UserInput.parameter_estimation_settings['multistart_initialPointsDistributionType'] = 'sobol'      
    UserInput.parameter_estimation_settings['multistart_parallel_sampling'] = True
    UserInput.parameter_estimation_settings['mcmc_exportLog'] = True #note that if we want the mcmc results for each parallel run to be exported, we need to state that, otherwise they won't be.
    
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)
    
    # run the model
    PE_object.doMultiStart()
    PE_object.createAllPlots() #This function calls each of the below functions so that the user does not have to.


