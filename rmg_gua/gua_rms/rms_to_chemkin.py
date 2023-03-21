

# %% [markdown]
# simple test script for executing rms from python, for use in the uncertainty pipeline. rms does sensitivities faster
# 

# %%
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
%matplotlib inline

# %%
file_dir = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/baseline/rms/chem53.rms"
phase_dict = rms.readinput(file_dir)

# %%
expt_condts = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/gua_cantera/all_experiments_reorg_sbr.yaml"

with open(expt_condts, 'r') as file:
    data = yaml.safe_load(file)

# pick just one experiment for example 
conditions = data[2]

# convert volume flow to molar flow
conditions["volume_flowrate"]

FC_temp = 293.15
conditions["molar_flow"] = conditions["volume_flowrate"] * 1.01325e5 / (8.3145 * FC_temp) 

# %%
# mechanism dictionaries index:  phase_dict[phasename]["Species" or "Reactions"]
gasspcs = phase_dict["gas"]["Species"]
gasrxns = phase_dict["gas"]["Reactions"]
surfacespcs = phase_dict["surface"]["Species"]
surfacerxns = phase_dict["surface"]["Reactions"]
interfacerxns = phase_dict[intfc_key]["Reactions"]

# %%
#Define the phase (how species thermodynamic and kinetic properties calculated)
ig = rms.IdealGas(gasspcs,gasrxns,name="gas") 
cat = rms.IdealSurface(surfacespcs, surfacerxns, 2.943e-5, name="surface")

# %%
# Set simulation gas Initial Temp and Pressure
initialcondsgas = {
        "T":conditions["temperature"],
        "P":conditions["pressure"],
        "CO":conditions["species"]["CO"],
        "CO2":conditions["species"]["CO2"],
        "H2":conditions["species"]["H2"],
} 
# Define the domain (encodes how system thermodynamic properties calculated)
domaingas,y0gas,pgas = rms.ConstantTPDomain(phase=ig,initialconds=initialcondsgas,sensitivity=True)

# %%
# Set simulation surf Initial Temp and Pressure
V = conditions["volume"]
A = conditions["catalyst_area"]
initialconds = {
        "T":conditions["temperature"],
        "A":conditions["catalyst_area"],
        "X":cat.sitedensity*A
} 
# Define the domain (encodes how system thermodynamic properties calculated)
domaincat,y0cat,pcat = rms.ConstantTAPhiDomain(phase=cat,initialconds=initialconds,sensitivity=True);

# %% [markdown]
# ## make reactor, inlet and outlet
# - makes an anonymous function x->42, is that velocity in? need to check if it is velocity or volume flowrate
# - also, I think the ```phi``` refers to chemical potential, but I should check, I think constantTPhi is just const T for our case. 

# %%
initialcondsinlet = {
        "T":conditions["temperature"],
        "P":conditions["pressure"],
        "CO":conditions["species"]["CO"],
        "CO2":conditions["species"]["CO2"],
        "H2":conditions["species"]["H2"],
    }

# construct reactor
inter,pinter = rms.ReactiveInternalInterfaceConstantTPhi(domaingas,domaincat,interfacerxns,initialcondsinlet["T"],A);

# make inlet and outlet
inletgas = rms.Inlet(domaingas,initialcondsinlet,Main.eval("x->"+str(conditions["molar_flow"])))
outletgas = rms.Outlet(domaingas,Main.eval("x->"+str(conditions["molar_flow"])))

# %%
# Define domains and interfaces
domains = (domaingas,domaincat)
interfaces = [inter,inletgas,outletgas]

# create a reactor for the system
react,y0,p = rms.Reactor(domains,(y0gas,y0cat),(0.0,100),interfaces,(pgas,pcat,pinter)) # Create the reactor object

# %%
run the simulation
t1 = time.time()
sol = de.solve(react.ode,de.CVODE_BDF(),abstol=1e-20,reltol=1e-8)
t2 = time.time()
print("elapsed time for sim: ", t2-t1)

# %%
ssys = rms.SystemSimulation(sol,domains,interfaces,p)

# %%
# load the chemkin file for the mechanism
chemkin_file = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/baseline/chemkin/chem_annotated-gas.inp"
chemkin_surf_file = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/baseline/chemkin/chem_annotated-surface.inp"
chemkin_dict = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/rmg_gua/baseline/chemkin/species_dictionary.txt"
model = ReactionModel()

model.species, model.reactions = load_chemkin_file(
    chemkin_file,
    chemkin_dict, 
    surface_path=chemkin_surf_file,
    )

# %%
sens_times = [1e-2, 3e-2, 100]
sens_rxn_dict = {}
max_len = 0
# get the ordering of reactions that are sensitive to CH3OH at different time scales. 
# our ordering is from highest sensitivity to low. so, start at one. if it is 
# found to be the first most sensitive at the first time value, then the 5th most sensitive at the 
# second time value, then the tenth at the third time value then the score it gets is 
# 1*5*10 = 50

for time in sens_times:
    sens_rxns, rxn_sens = rms.getrxntransitorysensitivities(ssys, "CH3OH", time)
    counter = 1
    for (rxn, sens) in zip(sens_rxns, rxn_sens):
        if rms.getrxnstr(rxn) in sens_rxn_dict.keys():
            old_sens = sens_rxn_dict[rms.getrxnstr(rxn)][0]
            sens_rxn_dict[rms.getrxnstr(rxn)][0] = counter*old_sens
            sens_rxn_dict[rms.getrxnstr(rxn)][1].append(sens)
            sens_rxn_dict[rms.getrxnstr(rxn)][2].append(counter)
        else:
            rxn_adj = rms.getrxnadjlist(rxn)
            reac_smiles = [Species().from_adjacency_list(reac) for reac in rxn_adj[0]]
            prod_smiles = [Species().from_adjacency_list(prod) for prod in rxn_adj[1]]
            sens_rxn_dict[rms.getrxnstr(rxn)] = [counter, [sens],[counter], reac_smiles, prod_smiles]

        counter +=1
    if len(sens_rxns) > max_len:
        max_len = len(sens_rxns)
    
# go through values in sens rxn dict that have fewer than len(sens_times) values 
# and multiply by max_len + 1. 
# penalizing for only being sensitive at one time point
for rxn in sens_rxn_dict.keys():
    if len(sens_rxn_dict[rxn][1]) < len(sens_times):
        sens_rxn_dict[rxn][0] *= (max_len + 1)
        for i in range(len(sens_times) - len(sens_rxn_dict[rxn][1])):
            sens_rxn_dict[rxn][1].append(0)
            sens_rxn_dict[rxn][2].append(max_len + 1)
            
# sort the dictionary by the score
sens_rxn_dict = {k: v for k, v in sorted(sens_rxn_dict.items(), key=lambda item: item[1][0])}

# print the dict: 
for rxn_str, entry in sens_rxn_dict.items():
    print(rxn_str, entry[0])

# %%
# now match up the rms reaction sensitivities with the actual rxn in chemkin
# match species does forward and reverse. if there are multiple matches, it will
# print a warning.
match_list = []
sens_cmkn_dict = {}
for rxn_str, entry in sens_rxn_dict.items():
    for rxn in model.reactions:
        counter = 0
        if rxn.matches_species(entry[3],entry[4]):
            print("match : ", rxn_str, counter)
            match_list.append(rxn_str)
            sens_cmkn_dict[rxn_str] = (entry[0], rxn)
        counter += 1
        if counter >=2: 
            print("more than 1 match found for ", rxn_str)
if len(match_list) == len(sens_rxn_dict.keys()):
    print("all matches found")
else: 
    for rxn_str, entry in sens_rxn_dict.items():
        if rxn_str not in match_list:
            print("no match found for ", rxn_str)

# %%
# get the specific node for the sensitive reaction. find out the highest populated parent node 

# load the kinetics database
RMG_db_folder= "/Users/blais.ch/Documents/_01_code/RMG_env_1/RMG-database/"
# Specify the path to the families
families_dir = os.path.join(RMG_db_folder,"input","kinetics","families")
if not os.path.exists(families_dir):
    raise OSError(f'Path to rules does not exist:\n{families_dir}')

# Specify the path to the libraries
kinetic_libraries_dir = os.path.join(RMG_db_folder,"input","kinetics","libraries","Surface")
if not os.path.exists(kinetic_libraries_dir):
    raise OSError(f'Path to kinetic libraries does not exist:\n{kinetic_libraries_dir}')

# do not load training, this will create more rules. in future, may want to 
# create those rules for perturbation, not sure. 
kinetics_families = ['all']
kinetics_database = KineticsDatabase()
kinetics_database.load_families(
        path=families_dir,
        families=kinetics_families,)

# %%
# breadcrumb cjb need to determine which sensitive reactions we should make rules for. 
# start with only one. add a test for most sensitive reaction, if a rule already exists, go for second most sensitive.

from rmgpy.data.kinetics.family import TemplateReaction
# test for reaction: get the nodes that would generate this rxn
rxn = model.reactions[30]
family = rxn.get_source()
template = rxn.template

# kinetics_database.families[family].get_labeled_reactants_and_products(rxn.reactants, rxn.products)
kinetics_database.families[family].add_atom_labels_for_reaction(rxn)
template = kinetics_database.families[family].get_reaction_template(rxn)
# kinetics_database.families[family].add_entry()
kinetics_database.families[family].has_rate_rule(template)
source_info = kinetics_database.families[family].extract_source_from_comments(rxn)
# new_rule = kinetics_database.families[family].get_kinetics_for_template(template)
parent = source_info[1][1]["rules"][0][0]
grp = rxn
name = 'O-C;VacantSite1;VacantSite2'
# kinetics_database.families[family].add_entry(parent, grp, name)
rank = source_info[1][1]["rules"][0][0].rank + 1


# %%

data = model.reactions[30].kinetics

if isinstance(data, StickingCoefficient):
    data = StickingCoefficientBEP(
        # todo: perhaps make a method StickingCoefficient.StickingCoefficientBEP
        #  analogous to Arrhenius.to_arrhenius_ep
        A=deepcopy(data.A),
        n=deepcopy(data.n),
        alpha=0,
        E0=deepcopy(data.Ea),
        Tmin=deepcopy(data.Tmin),
        Tmax=deepcopy(data.Tmax),
        coverage_dependence=deepcopy(data.coverage_dependence),
    )
elif isinstance(data, SurfaceArrhenius):
    data = SurfaceArrheniusBEP(
        # todo: perhaps make a method SurfaceArrhenius.toSurfaceArrheniusBEP
        #  analogous to Arrhenius.to_arrhenius_ep
        A=deepcopy(data.A),
        n=deepcopy(data.n),
        alpha=0,
        E0=deepcopy(data.Ea),
        Tmin=deepcopy(data.Tmin),
        Tmax=deepcopy(data.Tmax),
        coverage_dependence=deepcopy(data.coverage_dependence),
    )
data

# %%
from rmgpy.data.base import Entry
from rmgpy.reaction import Reaction

index = 3
data = data
new_entry = Entry(
    index=index,
    label=';'.join([g.label for g in template]),
    item=Reaction(reactants=[g.item for g in template], products=[]),
    data=data,
    rank=rank,
    short_desc="Rate rule generated for uncertainty",
    long_desc="Rate rule generated for uncertainty",
)
# new_entry.data.comment = "From training reaction {1} used for {0}".format(
#     ';'.join([g.label for g in template]), entry.index)

# new_entry.data.A.value_si /= entry.item.degeneracy
try:
    kinetics_database.families[family].rules.entries[new_entry.label].append(new_entry)
except KeyError:
    kinetics_database.families[family].rules.entries[new_entry.label] = [new_entry]
# index += 1

# %%
kinetics_database.families[family].rules.entries

# %%

kinetics_database.families[family].rules.save(os.path.join(families_dir, family, 'rules_' + 'test' + '.py'))

# %%
kinetics_database.families[family].groups.entries

# %%
for i in dir(model.reactions[30]):
    if "Reaction" in i:
        print(i)

# %%



