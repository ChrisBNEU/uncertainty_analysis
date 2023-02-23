import collections
import sys
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

################################################################
# collection of functions to make input files for PEUQSE
# 1. a RMS/chemkin lookup dictionary (they use diff naming conventions)
################################################################

def make_spc(spc):
    """
    make an RMG object from the rms object
    """
    if len(spc.adjlist) > 0:
        rmg_spc = Species().from_adjacency_list(spc.adjlist)
    else:
        rmg_spc = Species().from_smiles(spc.smiles)

    return rmg_spc


def initial_val_unc(kinetics, num):
    """ 
    return the initial value, an estimate of parameter uncertainty, and upper
     as a tuple
    A-factors: if not sticking coeffient, +/- 1, return LogA. else, +/- 0.5
    Ea: +/- 0.3 eV. 

    Limits: 
    A (upper):  - 1 for sticking coefficient, 1e25 for others
    A (lower):  - 0 for sticking coefficient, 0 for others
    E (upper):  - 400 kJ/mol
    E (lower):  - 0 kj/mol
    """
    # logA if it is not sticking coefficient
    if kinetics["A"] > 1:
        A = math.log10(kinetics["A"])
        A_unc = 1
        A_lb = 0
        A_ub = 25
        label_A = f"A_log_{num}"
    elif kinetics["A"] <= 1:
        A = kinetics["A"]
        A_unc = 1
        A_lb = 0
        A_ub = 1
        label_A = f"A_stick_{num}"

    E = kinetics["Ea"]
    E_unc = 30000000  # J/kmol for cantera
    E_lb = 0
    E_ub = 400000000  # J/kmol
    label_E = f"E_{num}"

    return (A, E), (A_unc, E_unc), (A_lb, E_lb), (A_ub, E_ub), (label_A, label_E)

def make_lookup_dict(base_path):
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

    rmg_2_ck_path = os.path.join(os.path.realpath(os.path.dirname(__file__)),"rmg_2_ck_dict.yaml") 
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

def make_ct_peuq_input(path, overwrite=True, sens_specs=["CH3OH"]):
    """
    generate the input required for the kinetic model parameters in peuqse.
    """
    # now load the sensitive reactions
    results_path = os.path.join(path, "rms", "results.csv")
    if os.path.exists(results_path) and not overwrite:
        res_df = pd.read_csv(results_path)
    else: 
        res_df = run_multi_reac(path, sens_specs)

    
    res_df = res_df[res_df["use_for_opt"] == True]

    scols = {}
    for spec in sens_spcs: 
        scols[spec] = res_df.filter(like=f"sens to {spec}").columns


    # there are a number of empty spots in the results (i.e. no sensitivity?)
    # pd.dropna
    highest_sens_dict = {}
    for names, cols in scols.items():

        res_df_T = res_df[cols].T
        sum_dict = {}
        for name, row in res_df_T.iterrows():
            total = row.abs().sum()
            sum_dict[name] = total

        # sort dict by value to get most sensitive
        highest = max(sum_dict, key=sum_dict.get)
        highest_sens_dict[names] = [highest, float(sum_dict[highest])]

    # make yaml for highest sensitivity reactions for each species
    highest_sens_path = os.path.join(os.path.realpath(
        os.path.dirname(__file__)), "highest_sens_rms.yaml")
    with open(highest_sens_path, "w") as f:
        yaml.dump(highest_sens_dict, f)

    # translate the rms reactions to chemkin reactions
    rms_2_ck = make_lookup_dict(path)
    sens_rxn_dict= {}
    for spec in sens_spcs:
        # get the reaction string        
        sens_rxn_str = highest_sens_dict[spec][0].split(" ")[0]

        # if there are duplicates skip
        if sens_rxn_str not in sens_rxn_dict.keys():
            sens_reac = sens_rxn_str.split("<=>")[0].split("+")
            sens_prod = sens_rxn_str.split("<=>")[1].split("+")
            sens_reac_ct = [rms_2_ck[name] for name in sens_reac]
            sens_prod_ct = [rms_2_ck[name] for name in sens_prod]

            sens_rxn_dict[sens_rxn_str] = [sens_reac_ct, sens_prod_ct]

    # create output dict for peuqse
    ct_peuqyml = {}
    ct_peuqyml["initial_values"] = []
    ct_peuqyml["labels"] = []
    ct_peuqyml["lower_bounds"] = []
    ct_peuqyml["uncertainties"] = []
    ct_peuqyml["upper_bounds"] = []

    # load the cantera mechanism
    cti_path = os.path.join(path, "cantera", "chem_annotated.cti")
    gas = ct.Solution(cti_path, "gas")
    surf = ct.Interface(cti_path, "surface1", [gas])

    # iterate over the sensitive reactions in the chemkin mechanism
    for sens_rxn, (sens_reac_ct, sens_prod_ct) in sens_rxn_dict.items():
        ct_reac = collections.Counter(sens_reac_ct)
        ct_prod = collections.Counter(sens_prod_ct)
        for num, rxn in enumerate(surf.reactions()):   
            reaclist = []
            prodlist = []
            for reac, stoic in rxn.reactants.items():
                if stoic > 1:
                    reaclist.extend([reac]*int(stoic))
                else:
                    reaclist.append(reac)

            for reac, stoic in rxn.products.items():
                if stoic > 1:
                    prodlist.extend([reac]*int(stoic))
                else:
                    prodlist.append(reac)

            reactants = collections.Counter(reaclist)
            products = collections.Counter(prodlist)

            fwd = reactants == ct_reac and products == ct_prod
            rev = reactants == ct_prod and products == ct_reac

            if rev or fwd:
                print(rxn)
                reac_dict = {}
                reac_dict["A"] = rxn.rate.pre_exponential_factor
                reac_dict["Ea"] = rxn.rate.activation_energy

                (A, E), (A_unc, E_unc), (A_lb, E_lb), (A_ub, E_ub), (label_A, label_E) =  initial_val_unc(reac_dict, num)
                
                ct_peuqyml["initial_values"].extend([A, E])
                ct_peuqyml["labels"].extend([label_A, label_E])
                ct_peuqyml["lower_bounds"].extend([A_lb, E_lb])
                ct_peuqyml["uncertainties"].extend([A_unc, E_unc])
                ct_peuqyml["upper_bounds"].extend([A_ub, E_ub])
                break

    # write the output yaml
    ct_param_path = os.path.join(os.path.realpath(os.path.dirname(__file__)), "ct_initial_small.yaml")
    with open(ct_param_path, "w") as f:
        yaml.dump(ct_peuqyml, f)


# make_ct_expt_file

# run the script with the location of the base model specied
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("please specify path to base model")
        sys.exit(1)

    base_model = sys.argv[1]

    sens_spcs = ["CH3OH", "CC", "CH4"]
    make_ct_peuq_input(base_model, overwrite=True, sens_specs=sens_spcs)
    # make_ct_expt_file(base_model, overwrite=False)




