import sys
import os
import yaml
import re
import math
import numpy as np
from copy import deepcopy
from rmgpy.chemkin import load_chemkin_file
from rmgpy.rmg.model import ReactionModel
from rmgpy.species import Species
from rmgpy.kinetics import StickingCoefficientBEP, StickingCoefficient, SurfaceArrheniusBEP, SurfaceArrhenius
# for getting type of reaction
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction 
from rmgpy.data.kinetics.database import KineticsDatabase
prefix = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
print("uncertainty repo path: ", prefix)

# use yaml safe load to preserve order (python 3.6+)
def make_rmg_reac_config(rmg_path, results_path=False, check_ranges=True):
    # load rmg database
    print("loading rmg database")
    kdb_path = os.path.join(rmg_path, "RMG-database", "input", "kinetics")
    kdb = KineticsDatabase()
    kdb.load(kdb_path, families='surface', depositories=False)

    # first get the values for A, E0, and Alpha from all of the surface families
    rule_dict = {}
    rule_unc_dict = {}
    rule_lb_dict = {}
    rule_ub_dict = {}
    rule_guess_dict = {}

    # gather the values, uncertainties, and bounds for each rule
    for fname, family in kdb.families.items(): 
        for rule_name, rule in family.rules.entries.items():
            entry = rule[0].data
            # check for surface arrhenius vs sticking coefficient
            if entry.A.value_si > 1.0: 
                A_val = math.log10(entry.A.value_si)
                A_unc = 1
                A_lb = 0
                A_ub = 25
                A_guess = A_val

            # changing to be the same as the surface arrhenius
            # else:
            #     A_val = entry.A.value_si
            #     A_unc = -1
            #     A_lb = 0
            #     A_ub = 1
            #     A_guess = 0.5
            else:
                A_val = math.log10(entry.A.value_si)
                A_unc = 1
                A_lb = -6
                A_ub = 1
                A_guess = A_val
                
            E0_val = entry.E0.value_si
            E0_unc = 30000  # J/mol, convert in cantera to j/kmol
            E0_lb = 0
            E0_ub = 400000  # J/mol
            
            # peuqse will take a "-1" input and just use range if it falls outside of uncertainty
            if E0_val - E0_unc <= E0_lb:
                
                # bump up our guess to the uncertainty 
                E0_guess = E0_unc + E0_lb
                E0_ub = E0_unc*2 + E0_lb
                # then change out uncertainty so it uses range
                E0_unc = -1

            else: 
                E0_guess = E0_val
            
            # we will treat all alphas the same, with allowable values spanning 0 to 1
            alpha_val = entry.alpha.value_si
            alpha_unc = -1
            alpha_lb = 0
            alpha_ub = 1
            alpha_guess = 0.5

            # if A, alpha, or E0 vals +/- unc goes out of bounds, then we will reset to 
            if check_ranges:
                if A_val + A_unc > A_ub or A_val - A_unc < A_lb:
                    print("A out of bounds for ", fname, " : ", rule_name, " : ", A_val, " +/- ", A_unc)
                    # A_val = A_guess
                if E0_val + E0_unc > E0_ub or E0_val - E0_unc < E0_lb:
                    print("E0 out of bounds for ", fname, " : ", rule_name, " : ", E0_val, " +/- ", E0_unc)
                    # E0_val = E0_guess
                if alpha_val + alpha_unc > alpha_ub or alpha_val - alpha_unc < alpha_lb:
                    print("alpha out of bounds for ", fname, " : ", rule_name, " : ", alpha_val, " +/- ", alpha_unc)
                    # alpha_val = alpha_guess

            rule_dict[fname + " : " +  rule_name] = {
                'A': A_val, 
                'E0': E0_val,
                'alpha': alpha_val,
                }  
            rule_unc_dict[fname + " : " +  rule_name] = {
                'A': A_unc, 
                'E0': E0_unc,
                'alpha': alpha_unc,
                }
            rule_lb_dict[fname + " : " +  rule_name] = {
                'A': A_lb, 
                'E0': E0_lb,
                'alpha': alpha_lb,
                }
            rule_ub_dict[fname + " : " +  rule_name] = {
                'A': A_ub, 
                'E0': E0_ub,
                'alpha': alpha_ub,
                }
            rule_guess_dict[fname + " : " +  rule_name] = {
                'A': A_guess, 
                'E0': E0_guess,
                'alpha': alpha_guess,
                }
            
    if results_path:  
        rule_config_file = os.path.join(results_path, "rule_config.yaml")
        with open(rule_config_file, 'w') as f:
            yaml.safe_dump(rule_dict, f, sort_keys=False)

        rule_unc_config_file = os.path.join(results_path, "rule_unc_config.yaml")
        with open(rule_unc_config_file, 'w') as f:
            yaml.safe_dump(rule_unc_dict, f, sort_keys=False)

        rule_lb_config_file = os.path.join(results_path, "rule_lb_config.yaml")
        with open(rule_lb_config_file, 'w') as f:
            yaml.safe_dump(rule_lb_dict, f, sort_keys=False)

        rule_ub_config_file = os.path.join(results_path, "rule_ub_config.yaml")
        with open(rule_ub_config_file, 'w') as f:
            yaml.safe_dump(rule_ub_dict, f, sort_keys=False)
            
        rule_guess_config_file = os.path.join(results_path, "rule_guess_config.yaml")
        with open(rule_guess_config_file, 'w') as f:
            yaml.safe_dump(rule_guess_dict, f, sort_keys=False)


def make_be_peuq_input(results_path=False):
    """
    append the binding energy parameters to the peuqse input files
    """

    be_values = {'C': 0., 'O': 0., 'N': 0., 'H': 0., 'vdw': 0.}
    be_unc = {'C': 3e7, 'O': 3e7, 'N': 3e7, 'H': 3e7, 'vdw': 2e7}
    be_lb = {'C': -1e8, 'O': -1e8, 'N': -1e8, 'H': -1e8, 'vdw': -1e8}
    be_ub = {'C': 1e8, 'O': 1e8, 'N': 1e8, 'H': 1e8, 'vdw': 1e8}
    be_guess = {'C': 0, 'O': 0, 'N': 0, 'H': 0, 'vdw': 0}
    if results_path:
        be_values_file = os.path.join(results_path,  "be_values.yaml")
        with open(be_values_file, 'w') as f:
            yaml.safe_dump(be_values, f, sort_keys=False)
        
        be_unc_file = os.path.join(results_path,  "be_unc.yaml")
        with open(be_unc_file, 'w') as f:
            yaml.safe_dump(be_unc, f, sort_keys=False)
        
        be_lb_file = os.path.join(results_path,  "be_lb.yaml")
        with open(be_lb_file, 'w') as f:
            yaml.safe_dump(be_lb, f, sort_keys=False)
        
        be_ub_file = os.path.join(results_path,  "be_ub.yaml")
        with open(be_ub_file, 'w') as f:
            yaml.safe_dump(be_ub, f, sort_keys=False)

        be_guess_file = os.path.join(results_path,  "be_guess.yaml")
        with open(be_guess_file, 'w') as f:
            yaml.safe_dump(be_guess, f, sort_keys=False)

def trim_rule_file(results_path, rules_used):
    """
    remove rules from rule config files that are not in the mechanism
    rules_used: dictionary of the rules used in the mechanism with a reaction or
    list of reactions as the values. 
    """

    rule_config_file = os.path.join(results_path, "rule_config.yaml")
    with open(rule_config_file, 'r') as f:
        rule_dict = yaml.safe_load(f)
    rule_unc_config_file = os.path.join(results_path, "rule_unc_config.yaml")
    with open(rule_unc_config_file, 'r') as f:
        rule_unc_dict = yaml.safe_load(f)
    rule_lb_config_file = os.path.join(results_path, "rule_lb_config.yaml") 
    with open(rule_lb_config_file, 'r') as f:
        rule_lb_dict = yaml.safe_load(f)
    rule_ub_config_file = os.path.join(results_path, "rule_ub_config.yaml")
    with open(rule_ub_config_file, 'r') as f:
        rule_ub_dict = yaml.safe_load(f)


    rule_list = list(rule_dict.keys())
    for rule in rule_list:
        if rule not in rules_used:
            print("removing rule: ", rule, " from config files")
            rule_dict.pop(rule)
            rule_unc_dict.pop(rule)
            rule_lb_dict.pop(rule)
            rule_ub_dict.pop(rule)
    
    # save all the yaml files 
    with open(rule_config_file, 'w') as f:
        yaml.safe_dump(rule_dict, f, sort_keys=False)
    with open(rule_unc_config_file, 'w') as f:
        yaml.safe_dump(rule_unc_dict, f, sort_keys=False)
    with open(rule_lb_config_file, 'w') as f:
        yaml.safe_dump(rule_lb_dict, f, sort_keys=False)
    with open(rule_ub_config_file, 'w') as f:
        yaml.safe_dump(rule_ub_dict, f, sort_keys=False)    

def make_ck_reac_config(results_path=False, trim_rules=False):
    """
    make config file detailing the rule associated with each rmg reaction
    trim_rule_file: if true, remove rules from rule config files that are not in 
    the mechanism
    """
    # now load the sensitive reactions
    path = os.path.join(prefix, "rmg_gua", "baseline")

    # load the chemkin file
    chemkin_file = os.path.join(path, "chemkin", "chem_annotated-gas.inp")
    chemkin_surf_file = os.path.join(path, "chemkin", "chem_annotated-surface.inp")
    chemkin_dict = os.path.join(path, "chemkin", "species_dictionary.txt")

    # load the chemkin file
    model = ReactionModel()
    model.species, model.reactions = load_chemkin_file(
        chemkin_file,
        dictionary_path=chemkin_dict,
        surface_path=chemkin_surf_file,
        use_chemkin_names=True,
    )
    ck_rule_dict = []
    rules_used = {}
    for rxn in model.reactions:
        if rxn.is_surface_reaction():
            if isinstance(rxn, TemplateReaction):
                comment = rxn.kinetics.comment
                if "Estimated using template" in comment: 
                    if "Average of" in comment:
                        template_str = re.search(r'Average of \[(.+?)\]', comment).group(1)
                    else:
                        template_str = re.search(r'Estimated using template \[(.+?)\]', comment).group(1)
                elif "Estimated using an average for rate rule" in comment:
                    template_str = re.search(r'Estimated using an average for rate rule \[(.+?)\]', comment).group(1)
                elif "Exact match found for rate rule" in comment:
                    template_str = re.search(r'Exact match found for rate rule \[(.+?)\]', comment).group(1)
                else:
                    print(f"no template found for {rxn.to_labeled_str(use_index=True)}")

                # get rxn path degeneracy
                if "Multiplied by reaction path degeneracy" in comment:
                    deg = float(re.search(r'Multiplied by reaction path degeneracy (.+?)\.', comment).group(1))
                else:
                    deg = 1.0
                
                if rxn.kinetics.A.value_si <=1.0:
                    A_val = rxn.kinetics.A.value_si
                    rtype = "stick"
                else: 
                    A_val = math.log10(rxn.kinetics.A.value_si)
                    rtype = "arrhenius"
                Ea_val = rxn.kinetics.Ea.value_si

                reac_dict = {
                    "A":A_val,
                    "Ea":Ea_val,
                }

                dhrxn = rxn.get_enthalpy_of_reaction(298)
                source_str = rxn.family + " : " + template_str
                reactants = [f"{spc.label}({spc.index})" for spc in rxn.reactants]
                products = [f"{spc.label}({spc.index})" for spc in rxn.products]

                # make a list 
                ck_rule_dict.append({rxn.to_labeled_str(use_index=True):{
                    "reactants": reactants, "products":products ,"source": source_str, 
                    "dHrxn": dhrxn, "reac_dict": reac_dict, "rtype":rtype, "deg":deg}})

                # add rule to rules_used dict. if it already exists, add the rxn to the list
                if source_str in rules_used.keys():
                    rules_used[source_str].append(rxn.to_labeled_str(use_index=True))
                else:
                    rules_used[source_str] = [rxn.to_labeled_str(use_index=True)]

    if results_path:
        ck_rule_dict_path = os.path.join(results_path, "ck_rule_dict.yaml")
        with open(ck_rule_dict_path, 'w') as f:
            yaml.dump(ck_rule_dict, f)

    if trim_rules:
        trim_rule_file(results_path, rules_used)

    return ck_rule_dict
 
def make_be_config(results_path=False, return_test_spec=False):
    """
    makes a dict containing each species and it's corresponding bond order

    returns:
        bond_orders (dict): dict of species and their bond orders
        test_spec (dict): dict of bond orders used for testing (cabon double bond, 
        carbon single bond, etc.)
    """
    # now load the sensitive reactions
    path = os.path.join(prefix, "rmg_gua", "baseline")

    # load the chemkin file
    chemkin_file = os.path.join(path, "chemkin", "chem_annotated-gas.inp")
    chemkin_surf_file = os.path.join(path, "chemkin", "chem_annotated-surface.inp")
    chemkin_dict = os.path.join(path, "chemkin", "species_dictionary.txt")

    # load the chemkin file
    model = ReactionModel()
    model.species, model.reactions = load_chemkin_file(
        chemkin_file,
        dictionary_path=chemkin_dict,
        surface_path=chemkin_surf_file,
        use_chemkin_names=True,
    )

    bond_orders = {}

    # We want a test species for each element w/ 2 different bond orders
    test_spec = {
        "vdw": "",
        "CS": "", 
        "CD": "",
        "CT": "",
        "OS": "", 
        "OD": "",
        "NS": "",
        "ND": "",
        "HS": "",
    }
    for spec in model.species: 
        if spec.contains_surface_site():

            adatoms = spec.molecule[0].get_adatoms()
            surface_sites = []
            isvdw = False
            for atom in spec.molecule[0].atoms:
                if atom.is_surface_site():
                    surface_sites.append(atom)
            normalized_bonds = {'C': 0., 'O': 0., 'N': 0., 'H': 0., 'vdw': 0.}
            max_bond_order = {'C': 4., 'O': 2., 'N': 3., 'H': 1.,} 
            

            for site in surface_sites:
                numbonds = len(site.bonds)
                if numbonds == 0 and len(spec.molecule[0].atoms) > 1:
                    isvdw = True
                    normalized_bonds['vdw'] = 1.
                    pass
                elif numbonds == 0 and len(spec.molecule[0].atoms) == 1:
                    print("empty surface site, skipping")
                    pass
                else:
                    assert len(site.bonds) == 1, "Each surface site can only be bonded to 1 atom"
                    bonded_atom = list(site.bonds.keys())[0]
                    bond = site.bonds[bonded_atom]
                    if bond.is_single():
                        bond_order = 1.
                    elif bond.is_double():
                        bond_order = 2.
                    elif bond.is_triple():
                        bond_order = 3.
                    elif bond.is_quadruple():
                        bond_order = 4.

                    normalized_bonds[bonded_atom.symbol] += bond_order / max_bond_order[bonded_atom.symbol]


            spec_ct_label = f"{spec.label}({spec.index})"
            bond_orders[spec_ct_label] = normalized_bonds
            print(f"{spec.label}({spec.index}) bond orders:{bond_orders[spec_ct_label]}")
            
            # test species
            if normalized_bonds['vdw'] == 1. and len(test_spec['vdw']) == 0:
                test_spec['vdw'] = spec_ct_label
            elif normalized_bonds['C'] == 0.25 and len(test_spec['CS']) == 0:
                test_spec['CS'] = spec_ct_label
            elif normalized_bonds['C'] == 0.5 and len(test_spec['CD']) == 0:
                test_spec['CD'] = spec_ct_label
            elif normalized_bonds['C'] == 0.75 and len(test_spec['CT']) == 0:
                test_spec['CT'] = spec_ct_label
            elif normalized_bonds['O'] == 0.5 and len(test_spec['OS']) == 0:
                test_spec['OS'] = spec_ct_label
            elif normalized_bonds['O'] == 1. and len(test_spec['OD']) == 0:
                test_spec['OD'] = spec_ct_label
            elif np.isclose(normalized_bonds['N'], 0.3) and len(test_spec['NS']) == 0:
                test_spec['NS'] = spec_ct_label
            elif normalized_bonds['N'] == 1. and len(test_spec['ND']) == 0:
                test_spec['ND'] = spec_ct_label
            elif normalized_bonds['H'] == 1. and len(test_spec['HS']) == 0:
                test_spec['HS'] = spec_ct_label

        else: 
            print(f"{spec.label}({spec.index}) is a gas phase species")
            pass
    # perturb each species BE via lsr
    if results_path:
        thermo_pert_file = os.path.join(results_path, "thermo_pert.yaml")
        with open(thermo_pert_file, 'w') as f:
            yaml.safe_dump(bond_orders, f, sort_keys=False)

        be_test_file = os.path.join(results_path, "be_test_config.yaml")
        with open(be_test_file, 'w') as f:
            yaml.safe_dump(test_spec, f, sort_keys=False)
    
    # return dict of test_cases if needed 
    if return_test_spec: 
        return bond_orders, test_spec
    else:
        return bond_orders
