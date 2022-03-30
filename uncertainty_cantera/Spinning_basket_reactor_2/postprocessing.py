import os
import pandas as pd
import numpy as np
import csv
import cantera as ct
import matplotlib.pyplot as plt

import seaborn
plt.style.use('seaborn-white')

import matplotlib.cm as cm
from rmgpy import chemkin

###############################################################################
# useful functions for plotting
###############################################################################

def plot_rates_grabow(df):
    '''
    This function returns a 2x2 set of charts like the ones above, 
    for comparison with the grabow rates and coverages. 
    df - dataframe to load
    labels - list of labels to use. first is the "x" values, then the y
    '''
    h2_mol = [0.5, 0.75, 0.8, 0.95]
    grabow_labels = [
        "Methanol Production",
        "Water-Gas Shift",
        "CO Hydrogenation",
        "CO2 Hydrogenation",
        "H2O Production",]
    
    h2_label = "y(H2)"
    co_co2_label = "CO2/(CO+CO2)"
    fig, ax = plt.subplots(2,2,figsize=(15,15))
    
    axes = [(0,0), (0,1), (1,0), (1,1)]
    
    color_dict = { 0:"k", 1:"y", 2:"g", 3:"b", 4:"r"}

    for coord,h2 in enumerate(h2_mol):
        
        for color,label in enumerate(grabow_labels):
            # "slot" where we place chart
            slot = axes[coord]
            df[np.isclose(df["y(H2)"], h2)].plot(x=co_co2_label,
                                                 y=label,
                                                 ax=ax[slot], 
                                                 color=color_dict[color],
                                                )
            ax[slot].set_title(f'mole frac H2 = {h2}')
            ax[slot].autoscale(enable=True, axis='y')
            ax[slot].set_ylabel("Turn over frequency ($s^{-1}$)")

def plot_rates_rmg(df):
    '''
    This function returns a 2x2 set of charts like the ones above, 
    for comparison with the grabow rates and coverages. 
    df - dataframe to load
    labels - list of labels to use. first is the "x" values, then the y
    '''
    
    h2_mol = [0.5, 0.75, 0.8, 0.95]

    #         "Methanol Production"  :  CH3O* + H* -> CH3OH* + *
    #         "Water-Gas Shift"      :  OH* + CO* -> COOH* + *
    #         "CO Hydrogenation"     :  CO* + H* -> HCO* + *
    #         "CO2 Hydrogenation"    :  CO2* + H* -> HCO2* + *
    #         "H2O Production"       :  #1 - #2
    #                                   #1:  H2O*+*->OH*+H*
    #                                   #2:  COOH* + OH* -> CO2* + H2O*

    rmg_labels = [
        "Methanol TOF ($s^{-1}$)",
        "WGS Reaction TOF ($s^{-1}$)",
        "CO Hydrogenation TOF ($s^{-1}$)",
        "CO2 Hydrogenation TOF ($s^{-1}$)",
        "H2O TOF ($s^{-1}$)",
    ]
        
    h2_label = "X_h2 initial"
    co_co2_label = "CO2/(CO2+CO)"
    color_dict = { 0:"k", 1:"y", 2:"g", 3:"b", 4:"r"}
    fig, ax = plt.subplots(2,2,figsize=(15,15))
    axes = [(0,0), (0,1), (1,0), (1,1)]

    for coord,h2 in enumerate(h2_mol):
        
        for color,label in enumerate(rmg_labels):
            # "slot" where we place chart
            slot = axes[coord]
            df[np.isclose(df[h2_label], h2)].plot(x=co_co2_label,
                                                 y=label,
                                                 ax =ax[slot], 
                                                 color=color_dict[color],
                                                )
            ax[slot].set_title(f'mole frac H2 = {h2}')
            ax[slot].autoscale(enable=True, axis='y')
            ax[slot].set_ylabel("Turn over frequency ($s^{-1}$)")

def plot_covs_grabow(df):
    '''
    This function returns a 2x2 set of charts like the ones above, 
    for comparison with the rmg rates and coverages. 
    df - dataframe to load
    '''
    
    h2_mol = [0.5, 0.75, 0.8, 0.95]
    grabow_labels = [
        "vacant",
        "HCOO",
        "CH3O",
        "H",
        "OH",
    ]
    
    h2_label = "y(H2)"
    co_co2_label = "CO2/(CO+CO2)"
    fig, ax = plt.subplots(2,2,figsize=(15,15))
    
    axes = [(0,0), (0,1), (1,0), (1,1)]
    
    color_dict = { 0:"k", 1:"y", 2:"g", 3:"b", 4:"r"}

    for coord,h2 in enumerate(h2_mol):
        
        for color,label in enumerate(grabow_labels):
            # "slot" where we place chart
            slot = axes[coord]
            df[np.isclose(df["y(H2)"], h2)].plot(x=co_co2_label,
                                                 y=label,
                                                 ax=ax[slot], 
                                                 color=color_dict[color],
                                                )
            ax[slot].set_title(f'mole frac H2 = {h2}')
            ax[slot].autoscale(enable=True, axis='y')
            ax[slot].set_ylabel("site fraction")

def plot_covs_rmg(df, labels, gas=False):
    '''
    This function returns a 2x2 set of charts like the ones above, 
    for comparison with the grabow rates and coverages. 
    df - dataframe to load
    labels - list of labels to use. first is the "x" values, then the y
    '''
    
    h2_mol = [0.5, 0.75, 0.8, 0.95]
    
    rmg_labels = labels
        
    h2_label = "X_h2 initial"
    co_co2_label = "CO2/(CO2+CO)"
    color_dict = { 0:"k", 1:"y", 2:"g", 3:"b", 4:"r", 5:"deeppink"}
    fig, ax = plt.subplots(2,2,figsize=(15,15))
    axes = [(0,0), (0,1), (1,0), (1,1)]

    for coord,h2 in enumerate(h2_mol):
        
        for color,label in enumerate(rmg_labels):
            # "slot" where we place chart
            slot = axes[coord]
            df[np.isclose(df[h2_label], h2)].plot(x=co_co2_label,
                                                 y=label,
                                                 ax =ax[slot], 
                                                 color=color_dict[color],
                                                )
            ax[slot].set_title(f'mole frac H2 = {h2}')
            ax[slot].autoscale(enable=True, axis='y')
            if gas: 
                ax[slot].set_ylabel("mole fraction")
            else: 
                ax[slot].set_ylabel("site fraction")

###############################################################################

# if len(sys.argv) < 2:
#     raise ValueError("Incorrect usage. Must pass the cantera model file as an argument to this analysis script")

# if not os.path.exists(sys.argv[1]):
#     raise OSError(f"Path to the cantera model file does not exist: {sys.argv[1]}")

# cti_file_path = sys.argv[1]
cti_file_path = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/run_0001/cantera/chem_annotated.cti"

rmg_model_folder = os.path.dirname(cti_file_path.replace("cantera/", ""))
csv_path = os.path.join(rmg_model_folder, "ct_analysis.csv")

# get grabow and model dictionary files
model_dict_file = rmg_model_folder + "/chemkin/species_dictionary.txt"
grabow_dict_file = "../../species_data/species_dictionary.txt"

# define model and grabow dictionary
model_dict = chemkin.load_species_dictionary(model_dict_file)
grabow_dict = chemkin.load_species_dictionary(grabow_dict_file)

# make a dictionary to "translate" the names from the grabow model to ours
# irrespective of the naming convention.
spc_trans = {}
for name, entry in model_dict.items(): 
    for g_name, g_entry in grabow_dict.items():
        if entry.is_isomorphic(g_entry):
            # remove (#) so it is neater
            g_name_new = g_name.split("(", 1)[0]
            spc_trans.update({g_name_new :name})


###############################################################################
# plotting
###############################################################################
 
 # construct dataframes
 rmg
