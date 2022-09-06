###############################################################################
# Sevy Harris 2021-09-02
# Chris Blais 2022-06-02
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
from rmgpy.data.surface import MetalDatabase
from rmgpy import constants
from tqdm import tqdm  # this is for the progress bar cause copying stuff takes a while
import numpy as np
import logging

# import stuff for easily writing config files
import pickle
import yaml

# needed for checking type
from rmgpy.kinetics import  StickingCoefficient, StickingCoefficientBEP

# future improvement: 
# ideally, load an input file, have this selectively go through and find what 
# should be perturbed, get atomic BE values initially and use those to generate 
# perturbation range
import platform
print(platform.python_version())

def get_range(orig, pert, type):
    """
    get the max and min value for a given input +/-pert value 
    orig - original value of parameter (a-factor, Ea, etc)
    pert - the maximum perturbation that can be applied (+/- 0.15 for 
           alpha perturbation, for example)
    type - whether it is a:
           stick - sticking coefficient
           A     - arrhenius
           Ea    - activation energy (for library or family E0 value)
           alpha - BEP alpha value 
           BE_J  - Binding energy value in J/mol
           BE_eV - Binding energy value in eV
    
    returns a tuple with the midpoint, min and max values for the range 
    """

    if type is "stick":
        max_allow = 1.0
        min_allow = 0.0
        
    # A is perturbed exponentially, so max and min are the exponent. 
    # original value is changes to exp
    elif type is "A":
        max_allow = 30.0
        min_allow = -30.0
        orig = np.log10(orig)

    elif type is "Ea":
        max_allow = 1.0e6
        min_allow = 0.0
    
    elif type is "alpha":
        max_allow = 1.0
        min_allow = 0.0

    elif type is "BE_J":
        max_allow = 9.6e5
        min_allow = -9.6e5

    elif type is "BE_eV":
        max_allow = 10
        min_allow = -10

    else: 
        logging.error(f"{type} type not specified when determining perturbation range")

    max = orig+pert
    min = orig-pert

    # check if above max or min allowable values
    if max > max_allow: 
        logging.warning(f"{type} maximum changed from {max} to {max_allow}")
        max = max_allow
    if min < min_allow: 
        logging.warning(f"{type} minimum changed from {min} to {min_allow}")
        min = min_allow

    mid = min + (max - min)/2

    if not np.isclose(mid, orig, 1e-3, 1e-3): 
        logging.warning(f"{type} value midpoint was changed from {orig} to {mid}")

    if max <= min: 
        logging.error(f"max {max} is less than that min {min}")

    return (mid, min, max)

def perturb_value(mid, pert, range_val):
    """
    orig  - original value for A, Ea, alpha, BE
    pert  - perturbation to be applied (0 to 1, from sobol sequence)
    range_val - min and max value that can be applied. 0 is min, 1, is max

    returns the perturbed value
    """
    sobolmin = 0
    sobolmax = 1.0
    sobolspan = sobolmax-sobolmin

    min = range_val[0]
    max = range_val[1]
    span = max - min 
    pert_val = (((pert - sobolmin) * span) / sobolspan) + min

    return pert_val

def generate_perturbed_files(
    RMG_db_folder, 
    output_path,
    perts,
    M
    ):
    """
    generates M perturbed files for global uncertainty analysis
    RMG_db_folder - path to rmg database
    output_path   - where the output files will be saved (such as the sobol map)
    M             - number of perturbed files to generate. each set will be labelled with
                    "XXXX" where XXXX is the perturbation from 1-M 

    creates:
    sobol_map.pickle       - a pickle file for the sobol map in the output_path directory
    sobol_range_map.pickle - the ranges of the perturbations applied to each parameter
    perturbed files        - for each kinetic family, kinetic library, thermo group, and 
                             thermo library, the values within are perturbed appropriately and given a 
                             perturbation number, e.g. surfaceThermoPt111_0001.py 
    """
    # Perturbation ranges for values (reference value +/- value below)
    DELTA_ALPHA_MAX_KN = perts["DELTA_ALPHA_MAX_KN"]
    DELTA_E0_MAX_J_MOL_KN = perts["DELTA_E0_MAX_J_MOL_KN"]
    DELTA_A_MAX_EXP_KN = perts["DELTA_A_MAX_EXP_KN"]
    DELTA_STICK_MAX_KN = perts["DELTA_STICK_MAX_KN"]
    DELTA_N_MAX_KN  = perts["DELTA_N_MAX_KN"] 

    # Perturbation ranges for unknown values (reference value +/- value below)
    DELTA_ALPHA_MAX_UNK = perts["DELTA_ALPHA_MAX_UNK"] 
    DELTA_E0_MAX_J_MOL_UNK = perts["DELTA_E0_MAX_J_MOL_UNK"] 
    DELTA_A_MAX_EXP_UNK = perts["DELTA_A_MAX_EXP_UNK"] 
    DELTA_STICK_MAX_UNK = perts["DELTA_STICK_MAX_UNK"]
    DELTA_N_MAX_UNK  = perts["DELTA_N_MAX_UNK"]

    # Thermo Perturbations
    # BEs that are in the metal database are in eV. vdw (in thermo db) is in
    # J/mol
    DELTA_E0_H = perts["DELTA_E0_H"]
    DELTA_E0_C = perts["DELTA_E0_C"]
    DELTA_E0_O = perts["DELTA_E0_O"]
    DELTA_E0_N = perts["DELTA_E0_N"]
    DELTA_E0_VDW = perts["DELTA_E0_VDW"]
    
    # load the things we want to perturb from a yaml
    thermo_libraries=perts["thermo_libraries"]
    thermo_groups_to_perturb = perts["thermo_groups"]
    lib_entries_to_perturb = perts["lib_entries_to_perturb"]
    kinetics_libraries = perts["kinetic_libraries"]
    kinetics_families_dict = perts["kinetics_families"]

    # make a list of just kinetics families for convenience
    kinetics_families = list(kinetics_families_dict.keys())

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
    # Specify the path to the families
    families_dir = os.path.join(RMG_db_folder,"input","kinetics","families")
    if not os.path.exists(families_dir):
        raise OSError(f'Path to rules does not exist:\n{families_dir}')

    # Specify the path to the libraries
    kinetic_libraries_dir = os.path.join(RMG_db_folder,"input","kinetics","libraries","Surface")
    if not os.path.exists(kinetic_libraries_dir):
        raise OSError(f'Path to kinetic libraries does not exist:\n{kinetic_libraries_dir}')

    # Specify the path to the thermo library
    os.path.join(RMG_db_folder, "input", "thermo")
    thermo_library_path = os.path.join(RMG_db_folder, "input", "thermo")
    if not os.path.exists(thermo_library_path):
        raise OSError(f'Path to rules does not exist:\n{thermo_library_path}')

    # Specify the path to the metal binding energy library
    metal_path = os.path.join(RMG_db_folder,"input","surface")
    if not os.path.exists(metal_path):
        raise OSError(f'Path to rules does not exist:\n{metal_path}')

    # Load the databases 
    kinetics_database = KineticsDatabase()
    kinetics_database.load_families(
        path=families_dir,
        families=kinetics_families,
    )
    kinetics_database.load_libraries(
        kinetic_libraries_dir,
        libraries=kinetics_libraries
    )

    metal_database = MetalDatabase()
    metal_database.load(metal_path)

    thermo_database = ThermoDatabase()
    thermo_database.load(
        thermo_library_path,
        libraries=thermo_libraries,
        depository=False,
        surface=True)

    # requires a reference database for the groups because deepcopy doesn't 
    # work for some reason. 
    thermo_database_ref = ThermoDatabase() 
    thermo_database_ref.load(
        thermo_library_path,
        libraries=thermo_libraries,
        depository=False,
        surface=True)
    
    # make a map of the perturbations to pickle later
    sobol_perturb_map = {}
    sobol_range_map = {}

    # get dimensions of sobol sequence 
    sobol_size = np.shape(x_sobol[:,0])

    # Kinetic family entry mapping
    # append the range as well. this info will be stored in the pickle file. 
    # for each sobol map entry: 
    #    [0] column index
    #    [1] perturbed value 
    #    [2] range of perturbations 
    if len(kinetics_families) > 0:
        for family_key in kinetics_database.families:
            family = kinetics_database.families[family_key]

            for entry_key in family.rules.entries.keys():
                # this is a good place for changes when I refactor, I 
                # currently have added a list to determine if the index is 
                # a node that should be perturbed more
                is_stick = isinstance(family.rules.entries[entry_key][0].data, StickingCoefficientBEP)
                entry_index = family.rules.entries[entry_key][0].index

                if entry_index in kinetics_families_dict[family_key]:
                    DELTA_ALPHA_MAX = DELTA_ALPHA_MAX_UNK
                    DELTA_E0_MAX_J_MOL = DELTA_E0_MAX_J_MOL_UNK
                    DELTA_A_MAX_EXP = DELTA_A_MAX_EXP_UNK
                    DELTA_STICK_MAX = DELTA_STICK_MAX_UNK 
                    DELTA_N_MAX = DELTA_N_MAX_UNK

                    E0_ref = 100000
                    alpha_ref = 0.5
                    # if it is a sticking coefficient reaction, set A to 0.5

                else: 
                    DELTA_ALPHA_MAX = DELTA_ALPHA_MAX_KN
                    DELTA_E0_MAX_J_MOL = DELTA_E0_MAX_J_MOL_KN
                    DELTA_A_MAX_EXP = DELTA_A_MAX_EXP_KN
                    DELTA_STICK_MAX = DELTA_STICK_MAX_KN
                    DELTA_N_MAX = DELTA_N_MAX_KN

                    E0_ref = family.rules.entries[entry_key][0].data.E0.value_si
                    alpha_ref = family.rules.entries[entry_key][0].data.alpha.value_si

                # label alpha column
                label = family_key + '/' + entry_key + '/alpha'
                # create array of nans for sobol map, so if one is unassigned we'll get an error
                sobol_map[label] = (sobol_col_index, np.empty((sobol_size))*np.nan)
                orig = alpha_ref 
                pert = DELTA_ALPHA_MAX
                type = "alpha"
                sobol_range_map[label] = get_range(orig, pert, type)
                sobol_col_index += 1

                # label pre-exponential column
                label = family_key + '/' + entry_key + '/A'
                sobol_map[label] = (sobol_col_index, np.empty((sobol_size))*np.nan)

                orig = family.rules.entries[entry_key][0].data.A.value_si
                
                if is_stick:
                    type = "stick"
                    pert = DELTA_STICK_MAX
                    orig = 0.5
                else:
                    type = "A"
                    pert = DELTA_A_MAX_EXP
                    orig = family.rules.entries[entry_key][0].data.A.value_si

                sobol_range_map[label] = get_range(orig, pert, type)
                sobol_col_index += 1

                # label activation energy column
                label = family_key + '/' + entry_key + '/E0'
                # create array of nans for sobol map, so if one is unassigned we'll get an error
                sobol_map[label] = (sobol_col_index, np.empty((sobol_size))*np.nan)
                orig = E0_ref
                pert = DELTA_E0_MAX_J_MOL
                type = "Ea"
                sobol_range_map[label] = get_range(orig, pert, type)
                sobol_col_index += 1

    all_library_entries = {}
    # Kinetic Library entry mapping
    if len(kinetics_libraries) > 0:
        # use known value perturbations for libraries
        DELTA_ALPHA_MAX = DELTA_ALPHA_MAX_KN
        DELTA_E0_MAX_J_MOL = DELTA_E0_MAX_J_MOL_KN
        DELTA_A_MAX_EXP = DELTA_A_MAX_EXP_KN
        DELTA_STICK_MAX = DELTA_STICK_MAX_KN
        DELTA_N_MAX = DELTA_N_MAX_KN
        
        for klib_key in kinetics_database.libraries:
            kinetics_lib = kinetics_database.libraries[klib_key]
            for klib_entry_key in kinetics_lib.entries.keys():
                kinetics_lib_entry = kinetics_lib.entries[klib_entry_key]
                is_stick = isinstance(kinetics_lib_entry.data, StickingCoefficient)

                if kinetics_lib_entry.label in lib_entries_to_perturb or 'all' in lib_entries_to_perturb:
                    
                    # label activation energy column
                    label = klib_key + '/' + str(klib_entry_key) + '/' + kinetics_lib_entry.label + '/Ea'
                    sobol_map[label] = (sobol_col_index, np.empty((sobol_size))*np.nan)
                    orig = kinetics_lib_entry.data.Ea.value_si
                    pert = DELTA_E0_MAX_J_MOL
                    type = "Ea"
                    sobol_range_map[label] = get_range(orig, pert, type)
                    sobol_col_index += 1

                    # label A-factor column
                    label = klib_key + '/' + str(klib_entry_key) + '/' + kinetics_lib_entry.label + '/A'
                    sobol_map[label] = (sobol_col_index, np.empty((sobol_size))*np.nan)
                    
                    if is_stick:
                        type = "stick"
                        pert = DELTA_STICK_MAX
                        orig = 0.5
                    else:
                        type = "A"
                        pert = DELTA_A_MAX_EXP
                        orig = kinetics_lib_entry.data.A.value_si

                    sobol_range_map[label] = get_range(orig, pert, type)
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
        "C", 
        "O", 
        "H",
        "N",
        "Vdw",
    ]

    # specify the metal you are perturbing
    metal = perts['metal']

    for label in scaling_groups:
        sobol_map[label] = (sobol_col_index, np.empty((sobol_size))*np.nan)
        sobol_col_index += 1

    metal_vals_orig = metal_database.libraries["surface"].entries[metal].binding_energies
    sobol_range_map["C"] = get_range(metal_vals_orig['C'].value, DELTA_E0_C, "BE_eV")
    sobol_range_map["O"] = get_range(metal_vals_orig['O'].value, DELTA_E0_O, "BE_eV")
    sobol_range_map["H"] = get_range(metal_vals_orig['H'].value, DELTA_E0_H, "BE_eV")
    sobol_range_map["N"] = get_range(metal_vals_orig['N'].value, DELTA_E0_N, "BE_eV")
    sobol_range_map["Vdw"] = get_range(0, DELTA_E0_VDW, "BE_J")

    # Perturb the values in the kinetics library
    if len(kinetics_libraries) > 0:
        print("Creating kinetics library files")
        for klib_key in kinetics_database.libraries:
            kinetics_lib = kinetics_database.libraries[klib_key]
            kinetics_lib_ref = copy.deepcopy(kinetics_lib)

            # if all entries selected, unload list of all entries to 
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
                            
                            # perturb the A factor 
                            sobol_key = klib_key + '/' + str(klib_entry_key) + '/' + kinetics_lib_entry.label + '/A'
                            sobol_col_index = sobol_map[sobol_key][0]
                            is_stick = isinstance(kinetics_lib_entry.data, StickingCoefficient)

                            mid = sobol_range_map[sobol_key][0]
                            pert = x_sobol[i, sobol_col_index] 
                            range_val = sobol_range_map[sobol_key][1:]

                            # we perturb A from 0<0.5<1 no matter what
                            if is_stick:
                                # for now, perturb all sticking coefficients from 0 to 1
                                A_perturbed = perturb_value(mid, pert, range_val)

                            else: 
                                # raise A to the perturbed power.
                                pert_exp = perturb_value(mid, pert, range_val)
                                A_perturbed = 10**pert_exp 
                            
                            kinetics_lib_entry.data.A.value_si = A_perturbed
                            # write perturbed value to sobol key
                            sobol_map[sobol_key][1][i] = A_perturbed


                            # perturb Ea
                            Ea_ref = kinetics_lib_ref.entries[klib_entry_key].data.Ea.value_si
                            sobol_key = klib_key + '/' + str(klib_entry_key) + '/' + kinetics_lib_entry.label + '/Ea'
                            sobol_col_index = sobol_map[sobol_key][0]

                            mid = sobol_range_map[sobol_key][0]
                            pert = x_sobol[i, sobol_col_index] 
                            range_val = sobol_range_map[sobol_key][1:]

                            Ea_perturbed = perturb_value(mid, pert, range_val)
                            kinetics_lib_entry.data.Ea.value_si = Ea_perturbed

                            # write perturbed value to sobol key
                            sobol_map[sobol_key][1][i] = Ea_perturbed

                kinetics_lib.save(os.path.join(kinetic_libraries_dir, klib_key, 'reactions_' + str(i).zfill(4) + '.py'))

    # Perturb the values in the kinetics families
    # check if index is listed to use a 
    if len(kinetics_families) > 0:
        print("Generating kinetics family files")
        for family_key in kinetics_database.families:
            print(family_key)
            family = kinetics_database.families[family_key]
            family_ref = copy.deepcopy(family)
            for i in tqdm(range(M)):
                for entry_key in family.rules.entries.keys():
                    entry = family.rules.entries[entry_key]
                    
                    is_stick = isinstance(family.rules.entries[entry_key][0].data, StickingCoefficientBEP)

                    # this is a good place for changes when I refactor, I 
                    # currently have added a list to determine if the index is 
                    # a node that should be perturbed more

                    ###########################################################
                    # perturb A-factor
                    ###########################################################
                    sobol_key = family_key + '/' + entry_key + '/A'
                    sobol_col_index = sobol_map[sobol_key][0]

                    mid = sobol_range_map[sobol_key][0]
                    pert = x_sobol[i, sobol_col_index] 
                    range_val = sobol_range_map[sobol_key][1:]

                    # we perturb A from 0<0.5<1 no matter what
                    if is_stick:
                        # for now, perturb all sticking coefficients from 0 to 1
                        A_perturbed = perturb_value(mid, pert, range_val)

                    else: 
                        # raise A to the perturbed power.
                        pert_exp = perturb_value(mid, pert, range_val)
                        A_perturbed = 10**pert_exp 
                            

                    entry[0].data.A.value = A_perturbed

                    # write perturbed value to sobol key
                    sobol_map[sobol_key][1][i] = A_perturbed  

                    ###########################################################
                    # Perturb the alpha value
                    ###########################################################
                    sobol_key = family_key + '/' + entry_key + '/alpha'
                    sobol_col_index = sobol_map[sobol_key][0]

                    mid = sobol_range_map[sobol_key][0]
                    pert = x_sobol[i, sobol_col_index] 
                    range_val = sobol_range_map[sobol_key][1:]
                    
                    alpha_perturbed = perturb_value(mid, pert, range_val)
                    entry[0].data.alpha.value = alpha_perturbed

                    # write perturbed value to sobol key
                    sobol_map[sobol_key][1][i] = alpha_perturbed

                    ###########################################################
                    # Perturb the E0 value
                    ###########################################################
                    sobol_key = family_key + '/' + entry_key + '/E0'
                    sobol_col_index = sobol_map[sobol_key][0]

                    mid = sobol_range_map[sobol_key][0]
                    pert = x_sobol[i, sobol_col_index] 
                    range_val = sobol_range_map[sobol_key][1:]

                    E0_perturbed = perturb_value(mid, pert, range_val)
                    entry[0].data.E0.value_si = E0_perturbed

                    # write perturbed value to sobol key
                    sobol_map[sobol_key][1][i] = E0_perturbed
                    
                family.rules.save(os.path.join(families_dir, family_key, 'rules_' + str(i).zfill(4) + '.py'))

    # perturb the values in the metal database
    # for vdw, we need to do every value 
    if len(scaling_groups) > 0:
        scaling_atoms = [
            "C", 
            "O", 
            "H",
            "N",
        ]
        print("generating binding energy perturbation files")
        metal_vals = metal_database.libraries["surface"].entries[metal]
        metal_vals_ref = copy.deepcopy(metal_vals)
        for i in tqdm(range(0, M)):
            # Get perturbations and record
            delta_E0_dict = {}
            for atom in scaling_atoms:
                delta_E0_dict[atom] = perturb_value(
                    sobol_range_map[atom][0], 
                    x_sobol[i, sobol_map[atom][0]], 
                    sobol_range_map[atom][1:],
                    )
                metal_vals.binding_energies[atom].value = delta_E0_dict[atom]
                sobol_map[atom][1][i] = delta_E0_dict[atom]

            metal_database.libraries['surface'].save(os.path.join(metal_path, 'libraries', 'metal' + '_' + str(i).zfill(4)+ ".py"))

    # unfortunately we still need to perturb the vdw groups and libraries 
    # because we do not have scaling for vdw
    if len(thermo_libraries) > 0:
        print("generating thermo library files")
        for library_key in thermo_database.libraries:
            thermo_lib = thermo_database.libraries[library_key]
            thermo_lib_ref = copy.deepcopy(thermo_lib)
            for i in tqdm(range(0, M)):

                mid = sobol_range_map["Vdw"][0] 
                pert = x_sobol[i, sobol_map["Vdw"][0]]
                range_val = sobol_range_map["Vdw"][1:]

                delta_E0_vdw = perturb_value(mid, pert, range_val)
                sobol_map["Vdw"][1][i] = delta_E0_vdw

                for entry_key in thermo_lib.entries.keys():
                    delta_E0 = 0
                    entry = thermo_lib.entries[entry_key]
                    # Don't perturb the energy level if it's just a vacant site
                    if entry_key is 'vacant':
                        continue
                    if entry.item.is_isomorphic(Molecule().from_adjacency_list("1 X  u0 p0 c0")):
                        continue
                    
                    # check which atom it is bound through. only perturb vdw
                    site_list = entry.item.get_surface_sites()

                    # skip if it is not a vdw species
                    if any(len(list(site.bonds.keys())) > 0 for site in site_list):
                        continue
                    
                    delta_E0 = delta_E0_vdw

                    # Perturb the E0 value, which is a5 in the NASA polymial
                    if entry.data.poly1 is not None:
                        E0_ref = thermo_lib_ref.entries[entry_key].data.poly1.c5
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
                thermo_lib.save(os.path.join(thermo_library_path, 'libraries', library_key + '_' + str(i).zfill(4) + '.py'))

    # perturb the values in the Thermo groups
    if len(thermo_groups_to_perturb) > 0: 
        print("generating thermo group files")
        for i in tqdm(range(0, M)):
            # Get perturbations and record
            mid = sobol_range_map["Vdw"][0] 
            pert = x_sobol[i, sobol_map["Vdw"][0]]
            range_val = sobol_range_map["Vdw"][1:]

            delta_E0_vdw = perturb_value(mid, pert, range_val)

            for group_name in thermo_groups_to_perturb:
                for group_entry_name in thermo_database.groups[group_name].entries:
                    delta_E0 = 0
                    if not thermo_database.groups[group_name].entries[group_entry_name].data:
                        continue  # only perturb thermo group entries that have thermo data
                    if isinstance(thermo_database.groups[group_name].entries[group_entry_name].data, str):
                        continue
                    thermo_group_entry = thermo_database.groups[group_name].entries[group_entry_name]
                    thermo_group_entry_ref = thermo_database_ref.groups[group_name].entries[group_entry_name]
                    
                    # get surface sites for group entry
                    site_list = thermo_database.groups[group_name].entries[group_entry_name].item.get_surface_sites()
                    
                    # make sure we have a vdw species
                    if any(len(list(site.bonds.keys())) > 0 for site in site_list):
                        continue

                    delta_E0 = delta_E0_vdw
                    H298_ref = thermo_group_entry_ref.data.H298.value_si
                    H298_perturbed = H298_ref + delta_E0
                    thermo_group_entry.data.H298.value_si = H298_perturbed

                thermo_database.groups[group_name].save(os.path.join(thermo_library_path, 'groups', group_name + '_' + str(i).zfill(4) + '.py'))


    # save the sobol map: 
    sobol_map_file = open(output_path + "sobol_map.pickle", "wb")
    pickle.dump(sobol_map, sobol_map_file)
    sobol_map_file.close()


    # for global uncertainty, it is advantageous to have the ranges saved as 
    # well. 
    sobol_range_file = open(output_path + "sobol_range_map.pickle", "wb")
    pickle.dump(sobol_range_map, sobol_range_file)
    sobol_range_file.close()

    # I am not sure we need to do this, redundant since we load in a yaml for this. 
    # # for convenience, save the list of perturbed families as a yaml
    # perturb_dict = {
    #     "thermo_libraries" : thermo_libraries,
    #     "thermo_groups" : thermo_groups_to_perturb,
    #     "lib_entries_to_perturb" : lib_entries_to_perturb,
    #     "kinetic_libraries" : kinetics_libraries,
    #     "kinetics_families" : kinetics_families_dict,
    #     "metal" : metal,
    #     ""
    # }
    # perturb_yaml = output_path + "perturb_groups.yaml"
    # with open(perturb_yaml, 'w') as f:
    #     # use safe_load instead of load
    #     yaml.safe_dump(perturb_dict, f)

if __name__ == "__main__":

    # for testing
    unc_folder = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
    conda_path = unc_folder + "conda/"
    RMG_db_folder = unc_folder + "RMG-database/"
    output_path = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/uncertainty_output_folder/"
    M = 10
    pert_yaml = os.path.join(unc_folder, "rmg_gua", "example", "perturb_groups.yaml")

    # define the stuff we want to perturb per yaml file
    with open(pert_yaml, 'r') as file:
        perts = yaml.safe_load(file)

    generate_perturbed_files(
    RMG_db_folder, 
    output_path,
    perts,
    M,
    )