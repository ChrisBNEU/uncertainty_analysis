import os
from pyrms import rms
from julia import Main
import time 
import yaml
from copy import deepcopy
from rmgpy.chemkin import load_chemkin_file
from rmgpy.rmg.model import ReactionModel
from rmgpy.species import Species
from rmgpy.kinetics import StickingCoefficientBEP, StickingCoefficient, SurfaceArrheniusBEP, SurfaceArrhenius
from rmgpy.data.kinetics.database import KineticsDatabase

# creates a lookup table of the rms name and ck/cantera name. 
# in retrospect, I think the name is basically the species label attribute. but 
# since I am not 100% sure, I will keep this script.

def make_spc(spc):
    """
    make an RMG object from the rms object
    
    """
    if len(spc.adjlist) > 0:
        rmg_spc = Species().from_adjacency_list(spc.adjlist)
    else:
        rmg_spc = Species().from_smiles(spc.smiles)

    return rmg_spc



base_pth = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline"
# base_path = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/baseline"

# load rms model
rms_path = os.path.join(base_path, "rms", "chem53.rms")
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


dir(model.species[12])
model.species[12].to_cantera(use_chemkin_identifier=True).name



# now make the species for the chemkin mechanism. 
# match each species string with it's chemkin counterpart
rmg_2_ck_dict = {}
for species in model.species:
    for rms_name, spec in rms_spc_dict.items():
        if species.is_isomorphic(spec):
            ck_name = species.to_cantera(use_chemkin_identifier=True).name
            rmg_2_ck_dict[rms_name] = ck_name
with open("rmg_2_ck_dict.yaml", "w") as f:
    yaml.dump(rmg_2_ck_dict, f)





