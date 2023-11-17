import collections
import itertools
import sys
import numpy as np 
from rmg_gua.gua_rms.sbr import rms_sbr
import os
from pyrms import rms
from julia import Main
import time
import yaml
import pandas as pd
import math
import cantera as ct
from copy import deepcopy
from rmgpy.chemkin import load_chemkin_file
from rmgpy.rmg.model import ReactionModel
from rmgpy.species import Species
from rmgpy.kinetics import StickingCoefficientBEP, StickingCoefficient, SurfaceArrheniusBEP, SurfaceArrhenius
from rmgpy.data.kinetics.database import KineticsDatabase
from matplotlib import pyplot as plt

###############################################################################
# various utilities for working with our awful rms->cantera->rms conversions
###############################################################################

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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


def make_spc(spc):
    """
    make an RMG object from the rms object
    """
    if len(spc.adjlist) > 0:
        rmg_spc = Species().from_adjacency_list(spc.adjlist)
    else:
        rmg_spc = Species().from_smiles(spc.smiles)

    return rmg_spc

def make_lookup_dict(base_path, results_path=None):
    """
    make a dictionary to look up the chemkin name from the rms name
    """

    # load rms model
    rms_path = get_highest_rms_file(base_path)
    phase_dict = rms.readinput(rms_path)

    # load chemkin files
    ck_path = os.path.join(base_path, "chemkin")
    cmkn_path = os.path.join(ck_path, "chem_annotated-gas.inp")
    cmkn_surf_path = os.path.join(ck_path, "chem_annotated-surface.inp")
    cmkn_dict_path = os.path.join(ck_path, "species_dictionary.txt")

    # get species from rms model
    gas_spec = phase_dict["gas"]["Species"]
    surf_spec = phase_dict["surface"]["Species"]

    # load the chemkin file for the mechanism
    model = ReactionModel()

    model.species, model.reactions = load_chemkin_file(
        cmkn_path,
        cmkn_dict_path,
        surface_path=cmkn_surf_path,
    )

    # make a dict of rms species with name as key, rmg obj as value
    rms_spc_dict = {}
    for spc in gas_spec:
        rmg_spc = make_spc(spc)
        rms_spc_dict[spc.name] = rmg_spc
    for spe in surf_spec:
        rmg_spc = make_spc(spe)
        rms_spc_dict[spe.name] = rmg_spc

    # now make the species for the chemkin mechanism.
    # match each species string with it's chemkin counterpart
    rmg_2_ck_dict = {}
    for species in model.species:
        for rms_name, spec in rms_spc_dict.items():
            if species.is_isomorphic(spec):
                ck_name = species.to_cantera(use_chemkin_identifier=True).name
                rmg_2_ck_dict[rms_name] = ck_name

    if results_path:
        rmg_2_ck_path = os.path.join(results_path, "rmg_2_ck_dict.yaml")
    
    else:
        rmg_2_ck_path = os.path.join(os.path.realpath(
            os.path.dirname(__file__)), "rmg_2_ck_dict.yaml")
        
    with open(rmg_2_ck_path, "w") as f:
        yaml.dump(rmg_2_ck_dict, f)

    return rmg_2_ck_dict

def run_reactor(condts, path, sens_spcs):

    # initialize reactor
    sbr_ss = rms_sbr(
        file_path=path,
        reac_config=condts,
        rtol=1.0e-11,
        atol=1.0e-22,
    )

    results = sbr_ss.run_simulation(sens_spcs=sens_spcs)
    return results

def get_highest_rms_file(path):
    """
    get the rms model with the highest number 
    """
    # get highest numbered model in rms folder
    rms_path = os.path.join(path, "rms")
    rms_files = os.listdir(rms_path)
    rms_files = [f for f in rms_files if f.endswith(".rms")]
    highest_path = ""

    # get the highest numbered file
    for filerms in rms_files:
        if highest_path == "":
            highest_path = filerms
        num_path = int(filerms.split(".")[0].split("chem")[1])
        highest_num = int(highest_path.split(".")[0].split("chem")[1])
        if num_path > highest_num:
            highest_path = filerms

    full_rms_path = os.path.join(rms_path, highest_path)
    return full_rms_path


def run_multi_reac(path, sens_spcs):
    """
    wip make multiprocess loop to run rms reactor, return dataframe of results
    """

    rms_filepath = get_highest_rms_file(path)

    if os.path.exists("/work"):
        settings_path = "/work/blais.ch/meOH_repos/uncertainty_analysis/rmg_gua/gua_cantera/experiments_reorg_onlyopt.yaml"
    else:
        settings_path = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/gua_cantera/experiments_reorg_onlyopt.yaml"

    csv_path = os.path.join(path, "rms", "results.csv")

    # load the settings
    with open(settings_path, "r") as f:
        settings = yaml.load(f, Loader=yaml.FullLoader)

    result = []

    count = 0
    ttot1 = time.time()
    for condition in settings:
        print("running ", count)

        t1 = time.time()
        res = run_reactor(condition, rms_filepath, sens_spcs=sens_spcs)
        result.append(res)
        t2 = time.time()

        print(f"process took {t2-t1} seconds")

    df = pd.DataFrame(result)
    df.to_csv(csv_path)

    return df

def trim_expts(expt_data):
    """
    trim the experiment list to a subset we've determined to be important. 
    8 experiments with high/low T, P, H2/(co2+co) and co/co2 ratios
    """
    # flatten the species sub-dictionary
    flat_expt_data = []
    for i, expt in enumerate(expt_data):
        flat_output = {}
        for key, value in expt.items(): 
            if key == "species":
                    for keysub, valuesub in expt[key].items():
                        flat_output[key+" "+keysub] = valuesub
            else: 
                flat_output[key] = value
        flat_expt_data.append(flat_output)

    flat_expt_df = pd.DataFrame(flat_expt_data)

    # make new column for h2 to co+co2 ratio and co to co2
    flat_expt_df["h2/(co+co2)"] = flat_expt_df["species H2"]/(flat_expt_df["species CO"] + flat_expt_df["species CO2"])
    flat_expt_df["co/co2"] = flat_expt_df["species CO"]/flat_expt_df["species CO2"]

    # get min/max values for each parameter
    pminmax = flat_expt_df["pressure"].min(), flat_expt_df["pressure"].max()
    tminmax = flat_expt_df["temperature"].min(), flat_expt_df["temperature"].max()
    ratioh2minmax = flat_expt_df["h2/(co+co2)"].min(), flat_expt_df["h2/(co+co2)"].max()
    ratiocominmax = flat_expt_df["co/co2"].min(), flat_expt_df["co/co2"].max()

    # get set of experiments that have combinations of low and high for each parameter
    expt_combos = {"p": pminmax, "t": tminmax, "ratioh2": ratioh2minmax, "ratioco": ratiocominmax}
    expt_combos = list(itertools.product(*expt_combos.values()))

    # Remove criteria for h2 ratio, makes it too restrictive, and without it we
    # end up having a variety of h2 rations anyway
    rtol = 1.0
    peuqse_expts = pd.DataFrame()
    for expset in expt_combos: 
        next_series = flat_expt_df[np.isclose(flat_expt_df["pressure"], expset[0], rtol=0.5) & \
            (np.isclose(flat_expt_df["temperature"], expset[1], atol=20)) & \
            (np.isclose(flat_expt_df["co/co2"], expset[3], rtol=0.5))]
        if next_series.empty:
            print("No experiments found for", expset)
        else: 
            # remove indices that are already in the dataframe
            next_series = next_series[~next_series.index.isin(peuqse_expts.index)]
            peuqse_expts = pd.concat([peuqse_expts,next_series.head(1)])

    # they didn't do experiments with [low h2 low co/co2] or [high h2 high co/co2] 
    # convert dataframe to dictionary
    # remove the h2/co+co2 and co/co2 columns we made
    peuqse_expts = peuqse_expts.drop(columns=["h2/(co+co2)", "co/co2"])
    peuqse_expts_dict = peuqse_expts.to_dict(orient="records")

    # now unflatten the species sub-dictionary
    new_peuqse_expts_dict = []
    for i, expt in enumerate(peuqse_expts_dict):
        flat_output = {}
        for key, value in expt.items(): 
            if key.startswith("species "):
                keysub = key.split(" ")[1]
                if "species" not in flat_output.keys():
                    flat_output["species"] = {}
                flat_output["species"][keysub] = value
            else: 
                flat_output[key] = value
        new_peuqse_expts_dict.append(flat_output)

    return new_peuqse_expts_dict

    

def make_ct_expt_file(results_path, use_peuq_expts=False):
    """
    make the ct experiment file for peuqse
    if use_peuq_expts is True, then use a smaller subset of experiments
    """
    
    settings_path = os.path.join(repo_dir, "rmg_gua", "gua_cantera", "experiments_reorg_onlyopt.yaml")
    # open the yaml file and get the list of experiments
    with open(settings_path, "r") as f:
        expt_list = yaml.load(f, Loader=yaml.FullLoader)

    if use_peuq_expts:
        expt_list = trim_expts(expt_list)
    
    # filter out the experiments we don't want
    new_expt_yaml = []
    for expt in expt_list:
        if "use_for_opt" in expt.keys():
            if expt["use_for_opt"] == True and expt["temperature"] < 518:
                new_expt_yaml.append(expt)

    # save the new yaml
    new_expt_file_orig = os.path.join(results_path, "expt_data_orig.yaml")
    with open(new_expt_file_orig, "w") as f:
        yaml.dump(new_expt_yaml, f)

    # load the experimental data yaml so we get a dict of lists for each parameter
    expt_data = repackage_yaml(new_expt_file_orig)

    # we need a new file for this
    new_expt_file = os.path.join(results_path, "expt_data.yaml")
    with open(new_expt_file, "w") as f:
        yaml.dump(expt_data, f)

    # make the uncertainties file
    unc_yaml = {}
    unc_yaml["catalyst_area"] = len(expt_data["catalyst_area"])*[7]
    unc_yaml["pressure"] = len(expt_data["pressure"])*[1e3]
    unc_yaml["species_CO"] = (np.array(expt_data["species_CO"])*0.01).tolist()
    unc_yaml["species_CO2"] = (np.array(expt_data["species_CO2"])*0.01).tolist()
    unc_yaml["species_H2"] = (np.array(expt_data["species_H2"])*0.01).tolist()
    unc_yaml["species_H2O"] = (np.array(expt_data["species_H2O"])*0.01).tolist()
    unc_yaml["species_out_CH3OH"] = (np.array(expt_data["species_out_CH3OH"])*0.01).tolist()
    unc_yaml["species_out_H2O"] = (np.array(expt_data["species_out_H2O"])*0.01).tolist()
    unc_yaml["species_out_H2"] = (np.array(expt_data["species_out_H2"])*0.01).tolist()
    unc_yaml["species_out_CO"] = (np.array(expt_data["species_out_CO"])*0.01).tolist()
    unc_yaml["species_out_CO2"] = (np.array(expt_data["species_out_CO2"])*0.01).tolist()
    unc_yaml["temperature"] = len(expt_data["temperature"])*[0.05]

    # save the new yaml
    new_unc_file = os.path.join(results_path, "expt_unc.yaml")
    with open(new_unc_file, "w") as f:
        yaml.dump(unc_yaml, f)
        

