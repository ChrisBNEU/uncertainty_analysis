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
sys.path.insert(0, repo_dir)
import rmg_gua.gua_peuqse.ct_simulation as ct_simulation
from rmg_gua.gua_peuqse.runtime_utilities import get_all_param_lists, make_exp_data_lists
project_path = os.path.dirname(os.path.abspath(__file__))


def setup_userinput(project_path, use_ranges=False, reduce_space=None):
    """
    sets up the common user inputs for a peuqse run.
    """
    
    ct_simulation.sim_init(project_path)
    results_path = os.path.join(project_path, "config")
    peuq_path = os.path.join(project_path, "peuqse")

    if not os.path.exists(peuq_path):
        os.mkdir(peuq_path)
    
    # get the peuqse parameters
    param_dict = get_all_param_lists(results_path=results_path, reduce_space=reduce_space)     

    
    data_path = os.path.join(repo_dir, "rmg_gua", "gua_cantera")
    x_data, y_data, y_unc = make_exp_data_lists(results_path)
    
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    y_unc = np.array(y_data)
    print(f"length is {len(x_data[0])} in main")
    
    UserInput.directories['graphs'] = os.path.join(peuq_path, "graphs/")
    UserInput.directories['logs_and_csvs'] = os.path.join(peuq_path, "logs_and_csvs/")
    UserInput.directories['pickles'] = os.path.join(peuq_path, "pickles/")
    
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
    UserInput.model['InputParameterInitialGuess'] = param_dict["guess_list"]

    if reduce_space is not None: 
        UserInput.model['reducedParameterSpace'] = param_dict["reduce_list"]
    
    #Can optionally change the initial guess to be different from prior means.
    if use_ranges:
        UserInput.model['InputParameterPriorValues_lowerBounds'] = param_dict["lower_list"]
        UserInput.model['InputParameterPriorValues_upperBounds'] = param_dict["upper_list"]
    
    UserInput.model['simulateByInputParametersOnlyFunction'] = ct_simulation.simulation_function_wrapper

    return UserInput
