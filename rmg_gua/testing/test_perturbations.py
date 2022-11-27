# ### black box testing 
# check output from RMG run in a chemkin file and ensure perturbs were applied correctly
import rmgpy.chemkin as chmkn
import random
from torch.quasirandom import SobolEngine
from rmgpy.data.kinetics.database import KineticsDatabase
import pickle
import os
import yaml

# needed for checking type
from rmgpy.kinetics import  StickingCoefficient, StickingCoefficientBEP

# perturbed model path 
pert_path = "/scratch/blais.ch/methanol_results_2022_07_19/"
# pert_path = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/"

with open(f'{pert_path}sobol_map.pickle', 'rb') as f:
    sobol_map = pickle.load(f)


# Perturbation ranges for known values (reference value +/- value below)
DELTA_ALPHA_MAX_KN = 0.15
DELTA_E0_MAX_J_MOL_KN = 30000  # 3 eV is about 30000 J/mol
DELTA_A_MAX_EXP_KN = 1
DELTA_STICK_MAX_KN = 0.5 # for sticking coefficient, perturb from 0-1 but by exponent.
DELTA_N_MAX_J_MOL_KN  = 1 

# Perturbation ranges for unknown values (reference value +/- value below)
DELTA_ALPHA_MAX_UNK = 0.5 #start value at 0.5 and perturb from 0-1
DELTA_E0_MAX_J_MOL_UNK = 100000  # 1 eV is about 100,000 J/mol. make sure eV is never 0.
DELTA_A_MAX_EXP_UNK = 4 # for sticking coefficient, perturb from 0-1 but by exponent. 
DELTA_STICK_MAX_UNK = 0.5 # for sticking coefficient
DELTA_N_MAX_J_MOL_UNK  = 1  

# Thermo Perturbation values
DELTA_E0_H = 0.3
DELTA_E0_C = 0.3
DELTA_E0_O = 0.3
DELTA_E0_N = 0.3
DELTA_E0_VDW = 20000


## sobol sequence testing
# Create the pseudo randoms
N = 10000
sobol = SobolEngine(dimension=300, scramble=True, seed=100)
x_sobol = sobol.draw(N)

# load original model
orig_path = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/"
orig_chmkn_path = os.path.join(orig_path, "chemkin")
orig_file_surf = os.path.join(orig_chmkn_path,"chem_annotated-surface.inp")
orig_file_gas = os.path.join(orig_chmkn_path,"chem_annotated-gas.inp")
orig_file_dict = os.path.join(orig_chmkn_path,"species_dictionary.txt")
orig_model = chmkn.load_chemkin_file(orig_file_surf, dictionary_path=orig_file_dict)
orig_model_gas = chmkn.load_chemkin_file(orig_file_gas, dictionary_path=orig_file_dict)

# load original kinetics database
db_path = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/RMG-database"
families_dir = os.path.join(db_path,"input","kinetics","families")

# Load the databases 
orig_kdb = KineticsDatabase()
orig_kdb.load_families(
    path=families_dir,
    families="all",
)

# get random models from 0-5000
model_dict = {}
db_dict = {}

seed_n = 4
test_run = False
for i in range(1,2):
    if test_run:
        model = i
        model_str = str(model).zfill(4)
    else:
        random.seed(i+seed_n)
        model = random.randint(9500, 10000)
        model_str = str(model).zfill(4)
        
        # scale the db#
        db = model - 9500
        db_str = str(db).zfill(4)
        
        print("model # ", model_str) 

    # load model
    chmkn_path = f"{pert_path}run_{model_str}/chemkin/"
    file_surf = chmkn_path + "chem_annotated-surface.inp"
    file_gas = chmkn_path +  "chem_annotated-gas.inp"
    file_dict = chmkn_path +  "species_dictionary.txt"
    chmkn_model = chmkn.load_chemkin_file(file_surf, dictionary_path=file_dict)
    
    # load database
    db_path = f"{pert_path}db_{db_str}/"
    families_dir = os.path.join(db_path,"input","kinetics","families")
    
    # Load the databases 
    pert_kdb = KineticsDatabase()
    pert_kdb.load_families(
        path=families_dir,
        families="all",
    )
    
    model_dict[model_str] = chmkn_model
    db_dict[model_str] = pert_kdb


# simple function to check thermo
def get_diff(spec, pert_spec):
    orig = spec.thermo.get_enthalpy(400)/9.6e4
    pert = pert_spec.thermo.get_enthalpy(400)/9.6e4
    diff = orig - pert
    return diff

################################################################
# get the rxn indices that correllate to the required test cases
# e.g. family - sticking coefficient - unknown rule
# dict of tested categories
test_categ = {}

# all lib entries use known perturbation
test_categ['lib_stick_index'] = -1
test_categ['lib_arrh_index'] = -1

# # known family rules
# test_categ['fam_stick_index_k'] = -1
# test_categ['fam_arrh_index_k'] = -1

# # unknown family rules
# test_categ['fam_stick_index_u'] = -1
# test_categ['fam_arrh_index_u'] = -1

# load perturb families file
example_path = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/example/"
pert_yaml_file = example_path + "perturb_groups.yaml"

with open(pert_yaml_file, "r") as f:
    pert_dict = yaml.safe_load(f)

# list of perturbed families with "unknown" rate rules (made up) 
unkn_fams = []
for fam, unks in pert_dict['kinetics_families'].items():
    if len(unks) >= 1: 
        unkn_fams.append(fam)
    
# # get all of the reactions to be tested
# for rxn_num, rxn in enumerate(orig_model[1]):
#     if hasattr(rxn, 'library'):
#         if isinstance(rxn.kinetics, StickingCoefficient) and test_categ['lib_stick_index'] == -1:
#             print("lib_stick_index")
#             display(rxn.to_labeled_str())
#             test_categ['lib_stick_index'] = rxn_num
#         elif test_categ['lib_arrh_index'] == -1: 
#             print("lib_arrh_rxn")
#             display(rxn.to_labeled_str())
#             test_categ['lib_arrh_index'] = rxn_num
#     else:
#         if isinstance(rxn.kinetics, StickingCoefficient):
#             if rxn.family in unkn_fams and test_categ['fam_stick_index_u'] == -1:
#                 print("fam_stick_rxn_u")
#                 display(rxn.to_labeled_str())
#                 test_categ['fam_stick_index_u'] = rxn_num
#             elif test_categ['fam_stick_index_k'] == -1: 
#                 print("fam_stick_rxn_kn")
#                 display(rxn.to_labeled_str())
#                 test_categ['fam_stick_index_k'] = rxn_num
                
#         else: 
#             if rxn.family in unkn_fams and test_categ['fam_arrh_index_u'] == -1:
#                 print("fam_arrh_rxn_u")
#                 display(rxn.to_labeled_str())
#                 test_categ['fam_arrh_index_u'] = rxn_num
#             elif test_categ['fam_arrh_index_k'] == -1: 
#                 print("fam_arrh_rxn_kn")
#                 display(rxn.to_labeled_str())
#                 test_categ['fam_arrh_index_k'] = rxn_num   
#####################################################################

test_categ = {
    'lib_stick_index': (0,"CPOX_Pt/Deutschmann2006_adjusted/1/H2 + Pt + Pt <=> HX + HX/"),
    'lib_arrh_index': (1, "CPOX_Pt/Deutschmann2006_adjusted/27/HOX + Pt <=> HX + OX/"),
    'fam_stick_index_k': (-1, 'Surface_Adsorption_Dissociative/C2H6;VacantSite1;VacantSite2/'),
    'fam_arrh_index_k': (-1,'Surface_Dissociation/Combined;VacantSite/'),
    'fam_stick_index_u': (-1,'Surface_EleyRideal_Addition_Multiple_Bond/Adsorbate1;Gas/'),
    'fam_arrh_index_u': (-1,'Surface_Dissociation_Double_vdW/AdsorbateVdW;VacantSite/'),
}
for cat, rxn_index in test_categ.items(): 
    if rxn_index == -1: 
        print(f"warning: test case {cat} not tested")
        
################################################################

for model_num, pert_model in model_dict.items(): 
    
    print(f"testing model {model_num}")
    
    
    E_0_c = sobol_map["C"][1][int(model_num)]
    E_0_o = sobol_map["O"][1][int(model_num)]
    E_0_h = sobol_map["H"][1][int(model_num)]
    E_0_vdw = sobol_map["Vdw"][1][int(model_num)]
    E_0_n = sobol_map["N"][1][int(model_num)]

    
    # check thermo perturbations
    for num, spec_num in enumerate(orig_model[0]):
        if spec_num.smiles == 'O=[Pt]':
            oxygen_spec = spec_num
            pert_oxygen_spec = pert_model[0][num]
            print("oxygen orig: ",oxygen_spec.thermo.get_enthalpy(298)/9.6e4)
            print("oxygen pert: ",pert_oxygen_spec.thermo.get_enthalpy(298)/9.6e4)
            print("oxygen sobol_pert: ", E_0_o)
            display(oxygen_spec)
        if spec_num.smiles == '[Pt]' and len(spec_num.molecule[0].atoms) > 1:
            hydrogen_spec = spec_num
            pert_hydrogen_spec = pert_model[0][num]
            print("hydrogen orig: ",hydrogen_spec.thermo.get_enthalpy(298)/9.6e4)
            print("hydrogen pert: ",pert_hydrogen_spec.thermo.get_enthalpy(298)/9.6e4)
            print("hydrogen sobol_pert: ", E_0_h)
            display(hydrogen_spec)
        if spec_num.smiles == 'O=C=[Pt]':
            carbon_spec = spec_num
            pert_carbon_spec = pert_model[0][num]
            print("carbon orig: ",carbon_spec.thermo.get_enthalpy(298)/9.6e4)
            print("carbon pert: ",pert_carbon_spec.thermo.get_enthalpy(298)/9.6e4)
            print("carbon sobol_pert: ", E_0_c)
            display(carbon_spec) 
            
    line = "#" * 50
            

    # check kinetic perturbation
    for case, data in test_categ.items(): 
        rxn_index = data[0]
        sobol_str = data[1]
        print(line)
        print(line)
        if hasattr(orig_model[1][rxn_index], 'library'):
            display(orig_model[1][rxn_index].to_labeled_str())

            for rxn_i, rxn in enumerate(pert_model[1]):
                if rxn.is_isomorphic(orig_model[1][rxn_index]):
                    pert_rxn_index = rxn_i
                

            # check that the A-factor is perturbed to the correct value. for sticking 
            # coefficient reactions this means anywhere from 0-1, independent of original 
            # value
            a_orig = orig_model[1][rxn_index].kinetics.A.value_si
            a_pert = pert_model[1][pert_rxn_index].kinetics.A.value_si
            sobol_pert = sobol_map[f"{sobol_str}A"][1][int(model_num)]
            print(case, rxn_index, pert_rxn_index)
            
            
            print(line)
            print("library A-factor stick")
            print("original: {0:e} \nPerturbed: {1:e} \nperturbation_value: {2:e}".format(a_orig, a_pert, sobol_pert))


            # check that the E0 is perturbed correctly
            Ea_orig = orig_model[1][rxn_index].kinetics.Ea.value_si
            Ea_pert = pert_model[1][pert_rxn_index].kinetics.Ea.value_si
            sobol_pert = sobol_map[f"{sobol_str}Ea"][1][int(model_num)]

            print(line)
            print("library Ea")
            print("original: {0:e} \nPerturbed: {1:e} \nperturbation_value: {2:e}".format(Ea_orig, Ea_pert, sobol_pert))

            
        else: 
            # this tests that the families reflect the sobol_map. 
            # family test: 
            rxn_index = data[0]
            sobol_key_str = data[1]
            family = sobol_key_str.split("/")[0]
            node = sobol_key_str.split("/")[1]

            print(line)
            print(f'testing family: {family}')
            print(line)

            # to get Ea and alpha, we have to use both values with the deltaH rxn

            # family Ea
            E0_orig = orig_kdb.families[family].rules.entries[node][0].data.E0.value_si
            E0_pert = db_dict[model_num].families[family].rules.entries[node][0].data.E0.value_si
            sobol_pert = sobol_map[sobol_key_str + "E0"][1][int(model_num)]
            print(line)
            print("family Ea")
            print("original: {0} \nPerturbed: {1} \nperturbation_value: {2}".format(E0_orig, E0_pert, sobol_pert))


            # family alpha
            alpha_orig = orig_kdb.families[family].rules.entries[node][0].data.alpha
            alpha_pert = db_dict[model_num].families[family].rules.entries[node][0].data.alpha
            sobol_pert = sobol_map[sobol_key_str + "alpha"][1][int(model_num)]

            print(line)
            print("family alpha")
            print("original: {0} \nPerturbed: {1} \nperturbation_value: {2}".format(alpha_orig, alpha_pert, sobol_pert))
            
            # family A
            A_orig = orig_kdb.families[family].rules.entries[node][0].data.A.value_si
            A_pert = db_dict[model_num].families[family].rules.entries[node][0].data.A.value_si
            sobol_pert = sobol_map[sobol_key_str + "A"][1][int(model_num)]


            print(line)
            print("family A")
            print("original: {0} \nPerturbed: {1} \nperturbation_value: {2}".format(A_orig, A_pert, sobol_pert))