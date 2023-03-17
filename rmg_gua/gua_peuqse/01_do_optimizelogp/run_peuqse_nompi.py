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
sys.path.insert(0, '../../../')

run_test = False#


import ct_simulation  # function for peuquse to optimize

if __name__ == "__main__":
    # Just a simple example. The user can also put the values in directly into the runfile or extract from a csv, for example.
    print("running job")
    # load the experimental data yaml so we get a dict of lists for each parameter

    with open("expt_data.yaml", "r") as f:
        expt_data = yaml.load(f, Loader=yaml.FullLoader)
        
    # load uncertainty yaml
    with open("expt_unc.yaml", "r") as f:
        expt_unc = yaml.load(f, Loader=yaml.FullLoader)
    
    # build the x-data array
    x_data = []
    x_data.append(expt_data["catalyst_area"])
    x_data.append(expt_data["pressure"])
    x_data.append(expt_data["species_CO"])
    x_data.append(expt_data["species_CO2"])
    x_data.append(expt_data["species_H2"])
    x_data.append(expt_data["temperature"])
    x_data.append(expt_data["volume_flowrate"])
    x_data = np.array(x_data)
    print(f"length is {len(x_data[0])} in main")
    # build x uncertainties array
    x_unc = []
    x_unc.append(expt_unc["catalyst_area"])
    x_unc.append(expt_unc["pressure"])
    x_unc.append(expt_unc["species_CO"])
    x_unc.append(expt_unc["species_CO2"])
    x_unc.append(expt_unc["species_H2"])
    x_unc.append(expt_unc["temperature"])
    x_unc.append(expt_unc["volume_flowrate"])
    x_unc = np.array(x_unc)


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

    #Provide the observed X values and Y values and uncertainties -- all should be arrays or lists with nesting like [[1,2,3]] or [[1,2,3,4],[4,5,6,6]]
    UserInput.responses['responses_abscissa'] = x_data
    UserInput.responses['responses_observed'] = y_data

    # estimated graaf uncertainties. 
    UserInput.responses['responses_observed_uncertainties'] = y_unc

    # load the input file we made for the imulation data (A, Ea, uncertainty A&Ea, etc.)
    ct_file = "./ct_initial_small.yaml" 
    with open(ct_file, "r") as f:
        mech_values = yaml.load(f, Loader=yaml.FullLoader)

    # convert to lists
    param_name_dict = {key: key for key in mech_values["labels"]} # dictionary of parameter names, not changing for now


    #Optional: provide labels for the responses axes and parameter names.
    UserInput.simulated_response_plot_settings['x_label'] = [
        "catalyst_area", "pressure (Pa)", "mole frac CO", "mole frac CO2", "mole frac H2", "temperature (K)", "volume flowrate (m^3/s)"]
    UserInput.simulated_response_plot_settings['y_label'] = [
        "CH3OH mole Frac out", "CO mole Frac out", "CO2 mole Frac out", "H2 mole Frac out", "H2O mole Frac out"]
    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = param_name_dict

    #Provide the prior distribution and uncertainties of the individual parameters.
    # prior expected values for a and b
    UserInput.model['InputParameterPriorValues'] = mech_values["initial_values"]
    # required. #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    UserInput.model['InputParametersPriorValuesUncertainties'] = mech_values["uncertainties"]
    # UserInput.model['InputParameterInitialGuess'] = #Can optionally change the initial guess to be different from prior means.
    UserInput.model['InputParameterPriorValues_lowerBounds'] = mech_values["lower_bounds"]
    UserInput.model['InputParameterPriorValues_upperBounds'] = mech_values["upper_bounds"]
    

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
