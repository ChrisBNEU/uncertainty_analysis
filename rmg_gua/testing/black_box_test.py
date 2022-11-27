# ### black box testing 
# check output from RMG run in a chemkin file and ensure perturbs were applied correctly
import rmgpy.chemkin as chmkn
import random
from torch.quasirandom import SobolEngine
import pickle
import os
import yaml

# perturbed model path 
pert_path = "/scratch/blais.ch/methanol_results_2022_07_19/"
code_path = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/"
pert_groups = f"{code_path}example/perturb_groups.yaml"

with open(f'{pert_path}sobol_map.pickle', 'rb') as f:
    sobol_map = pickle.load(f)

with open(pert_groups, 'r') as f:
    pert_yaml = yaml.safe_load(f)

# Perturbation ranges FOR values (reference value +/- value below)
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

# Thermo Perturbations
DELTA_E0_H = 30000
DELTA_E0_C = 30000
DELTA_E0_O = 30000
DELTA_E0_N = 30000
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

# get random models from 0-5000
model_dict = {}

for i in range(0,10):
    random.seed(i)
    model = random.randint(0, 5000)
    model_str = str(model).zfill(4)

    # load model
    chmkn_path = f"{pert_path}run_{model_str}/chemkin/"
    file_surf = chmkn_path + "chem_annotated-surface.inp"
    file_gas = chmkn_path +  "chem_annotated-gas.inp"
    file_dict = chmkn_path +  "species_dictionary.txt"
    chmkn_model = chmkn.load_chemkin_file(file_surf, dictionary_path=file_dict)

    model_dict[model_str] = chmkn_model


# simple function to check thermo
def get_diff(spec, pert_spec):
    orig = spec.thermo.get_enthalpy(400)/9.6e4
    pert = pert_spec.thermo.get_enthalpy(400)/9.6e4
    diff = orig - pert
    return diff

for model_num, pert_model in model_dict.items():
    
    # check thermo perturbations
    for num, spec_num in enumerate(orig_model[0]):
        if spec_num.smiles == 'O=[Pt]':
            oxygen_spec = spec_num
            pert_oxygen_spec = pert_model[0][num]
            print("oxygen orig: ",oxygen_spec.thermo.get_enthalpy(400)/9.6e4)
            print("oxygen pert: ",pert_oxygen_spec.thermo.get_enthalpy(400)/9.6e4)
            # display(oxygen_spec)
        if spec_num.smiles == '[Pt]' and len(spec_num.molecule[0].atoms) > 1:
            hydrogen_spec = spec_num
            pert_hydrogen_spec = pert_model[0][num]
            # display(hydrogen_spec)
        if spec_num.smiles == 'O=C=[Pt]':
            carbon_spec = spec_num
            pert_carbon_spec = pert_model[0][num]
            # display(carbon_spec) 
            
    E_0_c = float(DELTA_E0_C - 2.0 * x_sobol[int(model_num),sobol_map["C"][0]] * DELTA_E0_C)/9.6e4
    E_0_o = float(DELTA_E0_O - 2.0 * x_sobol[int(model_num),sobol_map["O"][0]] * DELTA_E0_H)/9.6e4
    E_0_h = float(DELTA_E0_H - 2.0 * x_sobol[int(model_num),sobol_map["H"][0]] * DELTA_E0_H)/9.6e4
    E_0_vdw = float(DELTA_E0_VDW  - 2.0 * x_sobol[int(model_num), sobol_map["Vdw"][0]] * DELTA_E0_VDW)/9.6e4
    E_0_n = float(DELTA_E0_N - 2.0 * x_sobol[int(model_num),sobol_map["N"][0]] * DELTA_E0_N)/9.6e4

    o_diff = get_diff(oxygen_spec, pert_oxygen_spec)
    h_diff = get_diff(hydrogen_spec, pert_hydrogen_spec)
    c_diff = get_diff(carbon_spec, pert_carbon_spec)

    print(model_num, " actual perturbation: ", o_diff, h_diff, c_diff)
    print(model_num, " intended perturb: ", E_0_o, E_0_h, E_0_c/2)

    # check kinetic perturbation
    # get first rxn from the mechanism
    rxn_index = 0
    # display(orig_model[1][rxn_index])

    # ensure it is the same as the first one in the perturbed model
    pert_model[1][0].is_isomorphic(orig_model[1][0])

    # check that the A-factor is perturbed to the correct value. for sticking 
    # coefficient reactions this means anywhere from 0-1, independent of original 
    # value
    a_orig = orig_model[1][rxn_index].kinetics.A
    a_pert = pert_model[1][rxn_index].kinetics.A
    sobol_pert = sobol_map["CPOX_Pt/Deutschmann2006_adjusted/1/H2 + Pt + Pt <=> HX + HX/A"][1][pert]

    print("library A-factor stick")
    print("original: {0} \nPerturbed: {1} \nperturbation_value: {2}".format(a_orig, a_pert, sobol_pert))


    # check that the E0 is perturbed correctly
    Ea_orig = orig_model[1][rxn_index].kinetics.Ea.value_si
    Ea_pert = pert_model[1][rxn_index].kinetics.Ea.value_si
    sobol_pert = sobol_map["CPOX_Pt/Deutschmann2006_adjusted/1/H2 + Pt + Pt <=> HX + HX/Ea"][1][pert]

    print("library Ea")
    print("original: {0} \nPerturbed: {1} \nperturbation_value: {2}".format(Ea_orig, Ea_pert, sobol_pert))

    # get first surface arrhenius rxn from the mechanism
    rxn_index = 7

    # ensure it is the same as the first one in the perturbed model
    pert_model[1][rxn_index].is_isomorphic(orig_model[1][rxn_index])
    print(type(pert_model[1][rxn_index].kinetics))
    print("rxn_source:", pert_model[1][rxn_index].family)

    sobol_col = 'CPOX_Pt/Deutschmann2006_adjusted/16/OCX + OX <=> CO2X + Pt/Ea'

    Ea_orig = orig_model[1][rxn_index].kinetics.Ea.value_si
    Ea_pert = pert_model[1][rxn_index].kinetics.Ea.value_si
    sobol_pert = sobol_map[sobol_col][1][pert]
    sobol_pert_value = x_sobol[pert,sobol_map[sobol_col][0]]
    delta_E0 = float(DELTA_E0_MAX_J_MOL_KN - 2.0 * sobol_pert_value * DELTA_E0_MAX_J_MOL_KN)

    print("library Ea stick rxn")
    print("original: {0} \nPerturbed: {1} \nperturbation_value: {2}".format(Ea_orig, Ea_pert, sobol_pert))
    print(delta_E0, Ea_pert-Ea_orig)


    sobol_col = 'CPOX_Pt/Deutschmann2006_adjusted/16/OCX + OX <=> CO2X + Pt/A'

    A_orig = orig_model[1][rxn_index].kinetics.A.value_si
    A_pert = pert_model[1][rxn_index].kinetics.A.value_si
    sobol_pert = sobol_map[sobol_col][1][pert]
    sobol_pert_value = x_sobol[pert,sobol_map[sobol_col][0]]
    delta_A = float(DELTA_A_MAX_EXP_KN - 2.0 * sobol_pert_value * DELTA_A_MAX_EXP_KN)
    A_perturbed = A_orig*10**delta_A

    print("library A surface arrhenius")
    print("original: {0:e} \nPerturbed: {1:e} \nperturbation_value: {2:e}".format(A_orig, A_pert, sobol_pert))
    print(A_perturbed, A_pert-A_orig)
    print(orig_model[1][rxn_index].kinetics.A)
    print(pert_model[1][rxn_index].kinetics.A)


    delta_A = float(DELTA_A_MAX_EXP_KN - 2.0 * sobol_pert_value * DELTA_A_MAX_EXP_KN)
    "{0:e}".format(A_orig*10**delta_A)

