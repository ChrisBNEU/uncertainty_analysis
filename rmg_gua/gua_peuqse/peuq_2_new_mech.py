# take the peuquse output and make it into a new mechanism
# will do both cantera and rms just to be complete. 
import pickle
from pyrms import rms
from diffeqpy import de
from julia import Main
import yaml
from julia import Sundials
from diffeqpy import de
import time
import matplotlib
from copy import deepcopy
from rmgpy.chemkin import load_chemkin_file
from rmgpy.rmg.model import ReactionModel
from rmgpy.species import Species
from rmgpy.kinetics import StickingCoefficientBEP, StickingCoefficient, SurfaceArrheniusBEP, SurfaceArrhenius
from rmgpy.data.kinetics.database import KineticsDatabase
import os
import copy
from peuqse_utilities import *



def get_peuq_map(path):
    """
    get the best fit parameters from peuqse
    """
    # path to peuqse pickle
    peuq_pickle = os.path.join(
        prefix, "rmg_gua", "gua_peuqse", "02_do_optimizelogp", 
        "pickles", "multistart_MAP_logP_and_parameters_values.pkl"
    )

    # load pickle 
    with open(peuq_pickle, "rb") as f:
        peuq_results = pickle.load(f)

    # sort the results array by the first value (logP)
    peuq_results = peuq_results[peuq_results[:,0].argsort()]

    # get the parameter labels from out input file
    ct_param_yaml = os.path.join(
        prefix, "rmg_gua", "gua_peuqse", "ct_initial_small.yaml"
    )
    with open(ct_param_yaml, "r") as f:
        ct_yaml = yaml.load(f, Loader=yaml.FullLoader)

    # make a dictionary of the parameters and their optimized values
    param_dict = {}
    for num, label in enumerate(ct_yaml["labels"]):
        # last row (-1) is max logp
        param_dict[label] = peuq_results[-1][i]

    # first, load the rms_file and change the parameters
    model_path = os.path.join(prefix, "rmg_gua", "baseline",)
    rms_file = get_highest_rms_file(model_path)
    phase_dict = rms.readinput(rms_file)

    # next, load the cantera objects
    ct_file = os.path.join(prefix, "rmg_gua", "baseline", "chem_annotated.cti")
    gas = ct.Solution(ct_file, "gas")
    surf = ct.Interface(ct_file, "surface1", [gas])
    ct_reactions = {}
    
    for label in param_dict.keys():
        num = int(label.split("_")[-1])
        eqn = surf.reactions()[num].equation

        # or param, value in new_rate_dict.items():
        new_rxn = self.surf.reactions()[num]

        A_i = self.surf.reactions()[num].rate.pre_exponential_factor
        Ea_i = surf.reactions()[num].rate.activation_energy
        b_i = surf.reactions()[num].rate.temperature_exponent

        if "A" in param and "log" in param:
            A_i = 10**float(value)
        elif "A" in param and "stick" in param:
            A_i = float(value)
        elif "E" in param:
            # input should use J/kmol
            Ea_i = float(value)
        else:
            logging.error(f"key {param} not recognized")

        rate = ct.Arrhenius(A=A_i, E=Ea_i, b=b_i)
        new_rxn.rate = rate
        print("new rxn rate", rate)
        surf.modify_reaction(num, new_rxn)

        print("newrxn: ", self.surf.reactions()[num].rate)
        if eqn not in ct_reactions.keys():
            ct_reactions[eqn] = surf.reactions()[num]

        # save as a pickle 

    # get the cantera reaction objects. match them to their rms reaction objects














if __name__ == "__main__":
    if os.path.exists("/work"):
        prefix = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
    else:
        prefix = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/"
    file_dir = os.path.join(prefix, "rmg_gua", "baseline", "rms", "chem32.rms")
    phase_dict = rms.readinput(file_dir)

    # save results to "00_run_fitted_model"
    updated_model_path = os.path.join("00_run_fitted_model", "updated_model")

