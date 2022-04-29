###############################################################################
# Sevy Harris 2021-09-02
# This script takes a handful of kinetics families, kinetics libraries, and
# thermo libraries and perturbs the E0, alpha, and Ea values
# Based off of Bjarne's Parametric Uncertainty Paper 10.1021/jacsau.1c00276
###############################################################################

import os
import copy
from torch.quasirandom import SobolEngine
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.molecule import Molecule
from rmgpy.data.thermo import ThermoDatabase
from rmgpy import constants
from tqdm import tqdm  # this is for the progress bar cause copying stuff takes a while
import numpy as np

# import stuff for easily writing config files
import pickle
import yaml

def generate_perturbed_files(
    RMG_db_folder, 
    output_path,
    M
    ):
    ###############################################################################
    # inputs
    ###############################################################################

    # directories to write to
    # database_path = "../RMG-database/"
    # results_path = "../uncertainty_output_folder/"

    # Perturbation ranges (reference value +/- value below)
    DELTA_ALPHA_MAX = 0.15
    DELTA_E0_MAX_J_MOL = 30000  # 3 eV is about 30000 J/mol
    DELTA_E0_MAX_J_MOL_VDW = 20000 # only perturb vdw species by +/- 0.2 eV

    # this will need to be appropriate units, and I have no idea what is default for RMG
    # a solution is just to put the max power that we want to raise the value to, 
    # then it will be perturbed by e4, e-4
    # after talking with Bjarne about this I think 1e1 is appropriate
    DELTA_A_MAX_EXP = 1

    # delta n? maybe -1 to 1 for now? this is basically
    # only here to fudge the numbers. might be useful for 
    # vdw  reactions (e.g. negative temperature dependence
    # on adsorption) 
    # DELTA_n_MAX_J_MOL = 1  

    # Define the number of perturbations to run
    # testing currently, use 10
    
    thermo_libraries = [
        'surfaceThermoPt111',
    ]

    # groups as in the group files (adsorptionpt111, ni111, etc.)
    thermo_groups_to_perturb = [
        'adsorptionPt111',
    ]

    # pick which entries to perturb in the kinetics library
    # WARNING: does not handle overlap of entries in different libraries
    # if 'all', perturb all entries
    lib_entries_to_perturb = [
        'all'
    ]
    # pick which kinetics libraries to perturb
    kinetics_libraries = [
        'CPOX_Pt/Deutschmann2006_adjusted'
    ]


    # Load the databases
    kinetics_families = [  # list the families to perturb
            "Surface_Abstraction",
            "Surface_Abstraction_Beta",
            "Surface_Abstraction_Beta_double_vdW",
            "Surface_Abstraction_Beta_vdW",
            "Surface_Abstraction_Single_vdW",
            "Surface_Abstraction_vdW",
            "Surface_Addition_Single_vdW",
            "Surface_Adsorption_Abstraction_vdW",
            "Surface_Adsorption_Bidentate",
            "Surface_Adsorption_Dissociative",
            "Surface_Adsorption_Dissociative_Double",
            "Surface_Adsorption_Double",
            "Surface_Adsorption_Single",
            "Surface_Adsorption_vdW",
            "Surface_Bidentate_Dissociation",
            "Surface_Dissociation",
            "Surface_Dissociation_Beta",
            "Surface_Dissociation_Beta_vdW",
            "Surface_Dissociation_Double",
            "Surface_Dissociation_Double_vdW",
            "Surface_Dissociation_vdW",
            "Surface_DoubleBond_to_Bidentate",
            "Surface_Dual_Adsorption_vdW",
            "Surface_EleyRideal_Addition_Multiple_Bond",
            "Surface_Migration",
            "Surface_vdW_to_Bidentate",
    ]

    ###############################################################################
    # Create Perturbed system
    ###############################################################################

    # Create the pseudo randoms
    sobol = SobolEngine(dimension=300, scramble=True, seed=100)
    x_sobol = sobol.draw(M)

    # make the map of sobol columns
    sobol_map = {}
    sobol_col_index = 0

    # make the map of sobol columns
    # make the map of sobol columns
    # Specify the path to the families
    families_dir = RMG_db_folder + "input/kinetics/families/"
    if not os.path.exists(families_dir):
        raise OSError(f'Path to rules does not exist:\n{families_dir}')

    # Specify the path to the libraries
    kinetic_libraries_dir = RMG_db_folder + "input/kinetics/libraries/Surface/"
    if not os.path.exists(kinetic_libraries_dir):
        raise OSError(f'Path to kinetic libraries does not exist:\n{kinetic_libraries_dir}')

    # Specify the path to the thermo library
    library_path = RMG_db_folder + "input/thermo/"
    if not os.path.exists(library_path):
        raise OSError(f'Path to rules does not exist:\n{library_path}')

    kinetics_database = KineticsDatabase()
    kinetics_database.load_families(
        path=families_dir,
        families=kinetics_families,
    )


    kinetics_database.load_libraries(
        kinetic_libraries_dir,
        libraries=kinetics_libraries
    )

    thermo_database = ThermoDatabase()
    thermo_database.load(
        library_path,
        libraries=thermo_libraries,
        depository=False,
        surface=True)

    thermo_database_ref = ThermoDatabase() # don't mess with this one
    thermo_database_ref.load(
        library_path,
        libraries=thermo_libraries,
        depository=False,
        surface=True)

    # make a map of the perturbations to pickle later
    sobol_perturb_map = {}

    # Kinetic family entry mapping
    if len(kinetics_families) > 0:
        for family_key in kinetics_database.families:
            family = kinetics_database.families[family_key]
            for entry_key in family.rules.entries.keys():
                # label BEP scaling parameter column
                label = family_key + '/' + entry_key + '/alpha'
                sobol_map[label] = (sobol_col_index,copy.deepcopy(x_sobol[:,sobol_col_index]))
                sobol_col_index += 1
                # label pre-exponential column
                label = family_key + '/' + entry_key + '/A'
                sobol_map[label] = (sobol_col_index,copy.deepcopy(x_sobol[:,sobol_col_index]))
                sobol_col_index += 1
                # label activation energy column
                label = family_key + '/' + entry_key + '/E0'
                sobol_map[label] = (sobol_col_index,copy.deepcopy(x_sobol[:,sobol_col_index]))
                sobol_col_index += 1

    all_library_entries = {}
    # Kinetic Library entry mapping
    if len(kinetics_libraries) > 0:
        for klib_key in kinetics_database.libraries:
            kinetics_lib = kinetics_database.libraries[klib_key]
            for klib_entry_key in kinetics_lib.entries.keys():
                kinetics_lib_entry = kinetics_lib.entries[klib_entry_key]
                if kinetics_lib_entry.label in lib_entries_to_perturb or 'all' in lib_entries_to_perturb:
                    
                    # label activation energy column
                    label = klib_key + '/' + str(klib_entry_key) + '/' + kinetics_lib_entry.label + '/Ea'
                    sobol_map[label] = (sobol_col_index,copy.deepcopy(x_sobol[:,sobol_col_index]))
                    sobol_col_index += 1

                    # label A-factor column
                    label = klib_key + '/' + str(klib_entry_key) + '/' + kinetics_lib_entry.label + '/A'
                    sobol_map[label] = (sobol_col_index,copy.deepcopy(x_sobol[:,sobol_col_index]))
                    sobol_col_index += 1

                    # if we have selected "all", then make a list of all 
                    # reactions so we know what to perturb later
                    if "all" in lib_entries_to_perturb:
                        if klib_key not in all_library_entries.keys():
                            all_library_entries[klib_key] = [kinetics_lib_entry.label]
                        else: 
                            all_library_entries[klib_key].append(kinetics_lib_entry.label)

    # in the case of correllated thermo uncertainties, we only have 4 parameters
    # dEh, dEo, dEc, and dEvdw, so we will have one column for each. 
    # added N_BE for the averaged nodes. but will likely not be an issue
    scaling_groups = [
        "H_BE", 
        "O_BE", 
        "C_BE",
        "N_BE",
        "vdw_BE",
    ]

    for label in scaling_groups:
        sobol_map[label] = (sobol_col_index,copy.deepcopy(x_sobol[:,sobol_col_index]))
        sobol_col_index += 1


    # Perturb the values in the kinetics library
    if len(kinetics_libraries) > 0:
        print("Creating kinetics library files")
        for klib_key in kinetics_database.libraries:
            kinetics_lib = kinetics_database.libraries[klib_key]
            kinetics_lib_ref = copy.deepcopy(kinetics_lib)

            #  if all entries selected, unload list of all entries to 
            # entries to perturb list
            if "all" in lib_entries_to_perturb:
                entries_to_perturb = all_library_entries[klib_key]
            else: 
                entries_to_perturb = lib_entries_to_perturb

            for i in tqdm(range(0, M)):
                for klib_entry_key in kinetics_lib.entries.keys():
                    kinetics_lib_entry = kinetics_lib.entries[klib_entry_key]

                    for label in entries_to_perturb:
                        if kinetics_lib_entry.label == label:

                            # perturb A-factor
                            # will need logic for sticking coefficient.
                            # for a=factor, stick coeff, perturb by a multiplier
                            # will need to account for units later
                            A_ref = kinetics_lib_ref.entries[klib_entry_key].data.A.value_si
                            sobol_key = klib_key + '/' + str(klib_entry_key) + '/' + kinetics_lib_entry.label + '/A'
                            sobol_col_index = sobol_map[sobol_key][0]
                            delta_A0 = (DELTA_A_MAX_EXP - 2.0 * x_sobol[i, sobol_col_index] * DELTA_A_MAX_EXP)
                            A_perturbed = A_ref*10**delta_A0
                            
                            kinetics_lib_entry.data.A.value_si = A_perturbed

                            # write perturbed value to sobol key
                            sobol_map[sobol_key][1][i] = delta_A0

                            # perturb Ea
                            Ea_ref = kinetics_lib_ref.entries[klib_entry_key].data.Ea.value_si
                            sobol_key = klib_key + '/' + str(klib_entry_key) + '/' + kinetics_lib_entry.label + '/Ea'
                            sobol_col_index = sobol_map[sobol_key][0]
                            delta_E0 = (DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[i, sobol_col_index] * DELTA_E0_MAX_J_MOL)
                            Ea_perturbed = Ea_ref + delta_E0
                            kinetics_lib_entry.data.Ea.value_si = Ea_perturbed

                            # write perturbed value to sobol key
                            sobol_map[sobol_key][1][i] = delta_E0

                kinetics_lib.save(os.path.join(kinetic_libraries_dir, klib_key, 'reactions_' + str(i).zfill(4) + '.py'))

    # Perturb the values in the kinetics families
    if len(kinetics_families) > 0:
        print("Generating kinetics family files")
        for family_key in kinetics_database.families:
            print(family_key)
            family = kinetics_database.families[family_key]
            family_ref = copy.deepcopy(family)
            for i in tqdm(range(M)):
                for entry_key in family.rules.entries.keys():
                    entry = family.rules.entries[entry_key]

                    # perturb A-factor
                    A_ref = family_ref.rules.entries[entry_key][0].data.A.value
                    sobol_key = family_key + '/' + entry_key + '/A'
                    sobol_col_index = sobol_map[sobol_key][0]
                    delta_A= DELTA_A_MAX_EXP - 2.0 * x_sobol[i, sobol_col_index] * DELTA_A_MAX_EXP

                    # raise A to the perturbed power.
                    A_perturbed = A_ref*10**delta_A
                    entry[0].data.A.value = A_perturbed

                    # write perturbed value to sobol key
                    sobol_map[sobol_key][1][i] = delta_A                

                    # Perturb the alpha value
                    alpha_ref = family_ref.rules.entries[entry_key][0].data.alpha.value
                    sobol_key = family_key + '/' + entry_key + '/alpha'
                    sobol_col_index = sobol_map[sobol_key][0]
                    delta_alpha = DELTA_ALPHA_MAX - 2.0 * x_sobol[i, sobol_col_index] * DELTA_ALPHA_MAX
                    alpha_perturbed = alpha_ref + delta_alpha
                    entry[0].data.alpha.value = alpha_perturbed

                    # write perturbed value to sobol key
                    sobol_map[sobol_key][1][i] = delta_alpha

                    # Perturb the E0 value
                    E0_ref = family_ref.rules.entries[entry_key][0].data.E0.value_si
                    sobol_key = family_key + '/' + entry_key + '/E0'
                    sobol_col_index = sobol_map[sobol_key][0]
                    delta_E0 = DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[i, sobol_col_index] * DELTA_E0_MAX_J_MOL
                    E0_perturbed = E0_ref + delta_E0
                    entry[0].data.E0.value_si = E0_perturbed

                    # write perturbed value to sobol key
                    sobol_map[sobol_key][1][i] = delta_E0
                    
                family.rules.save(os.path.join(families_dir, family_key, 'rules_' + str(i).zfill(4) + '.py'))

    # perturb the values in the Thermo library
    if len(thermo_libraries) > 0:
        print("generating thermo library files")
        for library_key in thermo_database.libraries:
            thermo_lib = thermo_database.libraries[library_key]
            thermo_lib_ref = copy.deepcopy(thermo_lib)
            for i in tqdm(range(0, M)):

                # Get perturbations and record
                delta_E0_C = DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[i, sobol_map["C_BE"][0]] * DELTA_E0_MAX_J_MOL
                delta_E0_O = DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[i, sobol_map["O_BE"][0]] * DELTA_E0_MAX_J_MOL
                delta_E0_H = DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[i, sobol_map["H_BE"][0]] * DELTA_E0_MAX_J_MOL
                delta_E0_N = DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[i, sobol_map["N_BE"][0]] * DELTA_E0_MAX_J_MOL
                delta_E0_vdw = DELTA_E0_MAX_J_MOL_VDW  - 2.0 * x_sobol[i, sobol_map["vdw_BE"][0]] * DELTA_E0_MAX_J_MOL_VDW 

                # replace sobol map value with the actual perturbation 
                sobol_map["C_BE"][1][i] = delta_E0_C
                sobol_map["O_BE"][1][i] = delta_E0_O
                sobol_map["H_BE"][1][i] = delta_E0_H
                sobol_map["N_BE"][1][i] = delta_E0_N
                sobol_map["vdw_BE"][1][i] = delta_E0_vdw

                for entry_key in thermo_lib.entries.keys():
                    delta_E0 = 0
                    entry = thermo_lib.entries[entry_key]
                    # Don't perturb the energy level if it's just a vacant site
                    if entry_key is 'vacant':
                        continue
                    if entry.item.is_isomorphic(Molecule().from_adjacency_list("1 X  u0 p0 c0")):
                        continue
                    
                    # check which atom it is bound through. bidentate
                    # perturbation will be additive (i.e. if it is bound 
                    # through O and C, it will add perturbation for each)
                    site_list = entry.item.get_surface_sites()
                    for site in site_list:
                        if len(list(site.bonds.keys())) > 0:
                            # currently we do not have one surface atom 
                            # bound to multiple adsorbate atoms
                            bonded_atom = list(site.bonds.keys())[0]
                            
                            if bonded_atom.is_carbon():
                                delta_E0 += delta_E0_C
                            elif bonded_atom.is_oxygen():
                                delta_E0 += delta_E0_O
                            elif bonded_atom.is_hydrogen():
                                delta_E0 += delta_E0_H
                            elif bonded_atom.is_nitrogen():
                                delta_E0 += delta_E0_N                    
                        # there are no instances where we have vdw and 
                        # a normal bond, so no +=
                        else: 
                            delta_E0 = delta_E0_vdw

                    # Perturb the E0 value, which is a5 in the NASA polymial
                    if entry.data.poly1 is not None:
                        E0_ref = thermo_lib_ref.entries[entry_key].data.poly1.c5
                        # E0_perturbed = E0_ref + delta_E0 / (constants.R / 1000.0)  # 8.314e-3
                        E0_perturbed = E0_ref + delta_E0 / constants.R  # 8.314
                        entry.data.poly1.c5 = E0_perturbed
                    if entry.data.poly2 is not None:
                        E0_ref = thermo_lib_ref.entries[entry_key].data.poly2.c5
                        E0_perturbed = E0_ref + delta_E0 / constants.R  # 8.314
                        entry.data.poly2.c5 = E0_perturbed
                    if entry.data.poly3 is not None:
                        E0_ref = thermo_lib_ref.entries[entry_key].data.poly3.c5
                        E0_perturbed = E0_ref + delta_E0 / constants.R  # 8.314
                        entry.data.poly3.c5 = E0_perturbed
                thermo_lib.save(os.path.join(library_path, 'libraries', library_key + '_' + str(i).zfill(4) + '.py'))

    # perturb the values in the Thermo groups
    if len(thermo_groups_to_perturb) > 0: 
        print("generating thermo group files")
        for i in tqdm(range(0, M)):
            # Get perturbations and record
            delta_E0_C = DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[i, sobol_map["C_BE"][0]] * DELTA_E0_MAX_J_MOL
            delta_E0_O = DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[i, sobol_map["O_BE"][0]] * DELTA_E0_MAX_J_MOL
            delta_E0_H = DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[i, sobol_map["H_BE"][0]] * DELTA_E0_MAX_J_MOL
            delta_E0_N = DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[i, sobol_map["N_BE"][0]] * DELTA_E0_MAX_J_MOL
            delta_E0_vdw = DELTA_E0_MAX_J_MOL_VDW  - 2.0 * x_sobol[i, sobol_map["vdw_BE"][0]] * DELTA_E0_MAX_J_MOL_VDW 

            for group_name in thermo_groups_to_perturb:

                for group_entry_name in thermo_database.groups[group_name].entries:
                    delta_E0 = 0
                    if not thermo_database.groups[group_name].entries[group_entry_name].data:
                        continue  # only perturb thermo group entries that have thermo data
                    if type(thermo_database.groups[group_name].entries[group_entry_name].data) is str:
                        continue
                    thermo_group_entry = thermo_database.groups[group_name].entries[group_entry_name]
                    thermo_group_entry_ref = thermo_database_ref.groups[group_name].entries[group_entry_name]
                    
                    # get surface sites for group entry
                    site_list = thermo_database.groups[group_name].entries[group_entry_name].item.get_surface_sites()

                    # check which atom it is bound through. bidentate
                    # perturbation will be additive (i.e. if it is bound 
                    # through O and C, it will add perturbation for each)
                    for site in site_list:
                        if len(list(site.bonds.keys())) > 0:
                            # currently we do not have one surface atom 
                            # bound to multiple adsorbate atoms
                            bonded_atom = list(site.bonds.keys())[0]
                            if bonded_atom.is_carbon():
                                delta_E0 += delta_E0_C
                            elif bonded_atom.is_oxygen():
                                delta_E0 += delta_E0_O
                            elif bonded_atom.is_nitrogen():
                                delta_E0 += delta_E0_N
                            elif bonded_atom.atomtype[0].label == "H": # no "is_hydrogen" method for groups :(
                                delta_E0 += delta_E0_H
                            elif "R" == bonded_atom.atomtype[0].label:
                                # average sobol perturbations for R*single_chemisorbed, since it is an averaged node. 
                                # this may need to be reevaluated for other wildcard nodes but currently 
                                # this is the only one. 
                                # averaged from C, O, and N nodes
                                sobol_avg = np.average([
                                    x_sobol[i, sobol_map["C_BE"][0]], 
                                    x_sobol[i, sobol_map["O_BE"][0]], 
                                    x_sobol[i, sobol_map["N_BE"][0]]]) # we are not doing nitrogen, use extra space in sobol
                                delta_E0+= DELTA_E0_MAX_J_MOL - 2.0 * sobol_avg * DELTA_E0_MAX_J_MOL

                            elif "R!H" == bonded_atom.atomtype[0].label:
                                # average sobol perturbations for R*bidentate, since it is an averaged node. 
                                # this may need to be reevaluated for other wildcard nodes but currently 
                                # this is the only one. 
                                # averaged from 3 C bidentate nodes and 1 N bidentate node. 
                                sobol_avg = np.average([
                                    x_sobol[i, sobol_map["C_BE"][0]], 
                                    x_sobol[i, sobol_map["C_BE"][0]], 
                                    x_sobol[i, sobol_map["C_BE"][0]],
                                    x_sobol[i, sobol_map["N_BE"][0]],]) # we are not doing nitrogen, use extra space in sobol
                                delta_E0+= DELTA_E0_MAX_J_MOL - 2.0 * sobol_avg * DELTA_E0_MAX_J_MOL
                        
                        # there are no instances where we have vdw and 
                        # a normal bond, so no +=
                        else: 
                            delta_E0 = delta_E0_vdw
                    
                    H298_ref = thermo_group_entry_ref.data.H298.value_si
                    H298_perturbed = H298_ref + delta_E0
                    thermo_group_entry.data.H298.value_si = H298_perturbed

                thermo_database.groups[group_name].save(os.path.join(library_path, 'groups', group_name + '_' + str(i).zfill(4) + '.py'))

    # save the sobol map: 
    sobol_map_file = open(output_path + "sobol_map.pickle", "wb")
    pickle.dump(sobol_map, sobol_map_file)
    sobol_map_file.close()

    # for convenience, save the list of perturbed families as a yaml
    perturb_dict = {
        "thermo_libraries" : thermo_libraries,
        "thermo_groups" : thermo_groups_to_perturb,
        "lib_entries_to_perturb" : lib_entries_to_perturb,
        "kinetic_libraries" : kinetics_libraries,
        "kinetics_families" : kinetics_families,
    }
    perturb_yaml = output_path + "perturb_groups.yaml"
    with open(perturb_yaml, 'w') as f:
        # use safe_load instead of load
        yaml.safe_dump(perturb_dict, f)
