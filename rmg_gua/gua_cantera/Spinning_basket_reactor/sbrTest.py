import cantera as ct
import math
import yaml
import numpy as np
import os
import shutil
import sys
from copy import deepcopy
import collections
from rmg_gua.gua_cantera.Spinning_basket_reactor.sbr import MinSBR
import pathlib
from rmg_gua.gua_cantera.Spinning_basket_reactor.make_peuq_config import make_ck_reac_config, make_rmg_reac_config, make_be_config

# tests for the sbr class 
class TestSBR:
    
    def setup_test(self, be_pert = False):
        """
        load an instance of the sbr class 
        load in data for rule and be perturbation
        load in config data for the reactions in the rule perturbation
        load in the species config for the be perturbation
        """    
        self.repo_dir = os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
        self.current_dir = os.path.join(
            self.repo_dir, "rmg_gua", "gua_cantera", "Spinning_basket_reactor") 

        # load the chemkin reaction config
        self.cti_path = os.path.join(
            self.repo_dir, "rmg_gua", "baseline", "cantera", "chem_annotated.cti")

        # if testing folder does not exist, create it
        if not os.path.exists(os.path.join(self.current_dir, "testing")):
            os.mkdir(os.path.join(self.current_dir, "testing"))
            
        self.results_path = os.path.join(self.current_dir, "testing")
        self.ck_rule_config = make_ck_reac_config(results_path = self.results_path)

        # load the rmg-database rule mapping
        self.rmg_path = os.path.dirname(os.environ["RMGPY"])

        # rmg database takes a while to load so if rmg rule file exists, 
        # just load the file. 
        if not os.path.exists(os.path.join(self.results_path, "rule_config.yaml")):
            print("re-load database since no file was found")
            self.rmg_rule_config = make_rmg_reac_config(self.rmg_path, self.results_path) 
        else: 
            print("rmg-database rule config file found, loading that")
            with open(os.path.join(self.results_path, "rule_config.yaml"), "r") as f:
                self.rmg_rule_config = yaml.load(f, Loader=yaml.FullLoader)

        # load the species binding energy mapping (which species has which element
        # bound to the surface)
        self.be_config, self.test_be_dict = make_be_config(
            results_path=self.results_path, return_test_spec=True)

        # make the perturbation list. If we are testing kinetics we do not want
        # to perturb the binding energy. 
        if be_pert:
            self.pert_list = self.make_test_pert_list(
                results_path=self.results_path, no_be_pert=False, rmg_rule_config=self.rmg_rule_config)
        else:
            self.pert_list = self.make_test_pert_list(
                results_path=self.results_path, no_be_pert=True, rmg_rule_config=self.rmg_rule_config)

        # create the reactor object
        conditions = {
            'catalyst_area': 44.424193000339805,
            'experiment_type': 'sbr',
            'expt_name': 'graaf_1988',
            'output': {'CH3OH': 0.00195844769830199, 'H2O': 0.0018249171734177633},
            'pressure': 3000000.0,
            'species': {'CO': 0.065, 'CO2': 0.261, 'H2': 0.674, 'H2O': 0},
            'species_out': {'CH3OH': 0.0044,
            'CO': 0.0662,
            'CO2': 0.2561,
            'H2': 0.6692,
            'H2O': 0.0041},
            'temperature': 483.5,
            'use_for_opt': True,
            'volume': 0.0001346957850226624,
            'volume_flowrate': 1.423e-05
            }
        
        self.sbr = MinSBR(
            self.cti_path,
            reac_config=conditions,
            rtol=1.0e-11,
            atol=1.0e-22,
            reaction_list=self.pert_list,
            results_path=self.results_path,
        ) 
        # create a copy of the reactor object with no perturbations
        self.orig_sbr = MinSBR(
            self.cti_path,
            reac_config=conditions,
            rtol=1.0e-11,
            atol=1.0e-22,
            reaction_list=None, 
            results_path=None,
        ) 
    def teardown(self):
        """
        remove the testing folder
        """
        shutil.rmtree(self.results_path)
        
    def make_test_pert_list(self, results_path=False, no_be_pert=False, rmg_rule_config=None):
        """
        make a test perturbation list 
        for testing, 2 lists are required. 
        1. perturbation list with different values for binding energy
        2. perturbation list with zero perturbation for binding energy, since 
        it affects the activation energy calculation. 
        """

        # load the perturbation list
        if no_be_pert:
            thermo_pert_dict = {'C': 0., 'O': 0., 'N': 0., 'H': 0., 'vdw': 0.}
        else:
            thermo_pert_dict = {'C': -3e7, 'O': -2e7, 'N': -1e7, 'H': 1e7, 'vdw': 2e7}

        self.thermo_pert_dict = thermo_pert_dict
        # add the rate rule perturbations 
        pert_list_test = []
        for rule in rmg_rule_config.values():
            pert_list_test.extend([rule["A"], rule["E0"], rule["alpha"]])

        # add the thermo perturbations
        pert_list_test.extend(list(thermo_pert_dict.values()))

        # write the perturbations to a file for additional troubleshooting
        if results_path:
            pert_list_path = os.path.join(results_path, "pert_list.yaml")
            with open(pert_list_path, 'w') as f:
                yaml.dump(pert_list_test, f)

        return pert_list_test

    def test_change_be(self):
        # test that perturbations are correct (should change by perturbation value at 298K)
        # does not acount for bidentates
        self.setup_test(be_pert=True) 

        kj_thermo_pert = {k: v/1e6 for k, v in self.thermo_pert_dict.items()}
        # go through the test dictionary for species matching the bond order 
        for bo, spec in self.test_be_dict.items():
            if len(spec) > 0:
                if bo == "vdw":
                    expected_diff = kj_thermo_pert["vdw"]
                elif bo == "CS":
                    expected_diff = kj_thermo_pert["C"] * 0.25
                elif bo == "CD":
                    expected_diff = kj_thermo_pert["C"] * 0.5
                elif bo == "CT":
                    expected_diff = kj_thermo_pert["C"] * 0.75
                elif bo == "OS":
                    expected_diff = kj_thermo_pert["O"] * 0.5
                elif bo == "OD":
                    expected_diff = kj_thermo_pert["O"] * 1.0
                elif bo == "NS":
                    expected_diff = kj_thermo_pert["N"] / 3
                elif bo == "ND":
                    expected_diff = kj_thermo_pert["N"] * (2/3)
                elif bo == "HS":
                    expected_diff = kj_thermo_pert["H"] * 1.0

                thermo_new = self.sbr.surf.species(self.sbr.surf.species_index(spec)).thermo
                thermo_old = self.orig_sbr.surf.species(self.orig_sbr.surf.species_index(spec)).thermo
                diff = thermo_new.h(298.15)/1e6 - thermo_old.h(298.15)/1e6
                print(f"compare {bo} pert to actual pert for {spec}")
                print("diff: ", diff)
                print("expected diff: ", expected_diff)
            
                assert np.isclose(diff, expected_diff, rtol=1e-3, atol=1e-3), f"{bo} perturbation is not correct"

    def test_change_reactions(self):
        """
        Test that the reaction rates are being changed correctly
        Nothing is enforced, but we at least print the output to see if it's reasonable
        """
        # setup class
        self.setup_test(be_pert=False) 

        ck_rule_dict_path = os.path.join(self.results_path,  "ck_rule_dict.yaml")
        with open(ck_rule_dict_path, 'r') as f:
            ck_rule_dict = yaml.load(f, Loader=yaml.FullLoader)

        rule_config_file = os.path.join(self.results_path, "rule_config.yaml")
        with open(rule_config_file, "r") as f:
            rule_dict = yaml.safe_load(f)
        ck_matches = {}
        ck_no_match = []
        for ck_num, ck_items in enumerate(ck_rule_dict):
            ck_rxn = list(ck_items.keys())[0] 
            ck_data = list(ck_items.values())[0]
            ck_reac = collections.Counter(ck_data["reactants"])
            ck_prod = collections.Counter(ck_data["products"])

            num = ck_num
            rxn = self.sbr.surf.reactions()[ck_num]
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

            fwd = reactants == ck_reac and products == ck_prod
            rev = reactants == ck_prod and products == ck_reac

            # if we have a match, check A, Ea against ck rxn
            if fwd:
                A_old = ck_data["reac_dict"]["A"]
                Ea_old = ck_data["reac_dict"]["Ea"]

                # record initial values from cantera
                A_i = self.sbr.surf.reactions()[num].rate.pre_exponential_factor
                Ea_i = self.sbr.surf.reactions()[num].rate.activation_energy

                n_reactants = len(ck_data["reactants"])

                if ck_data["rtype"] == "arrhenius":
                    A_old = 10**float(A_old)
                    # for bimolecular reactions we have m^2/mol*s in chemkin, 
                    # but m^2/kmol*s in cantera. all others are either sticking 
                    # (trimolecular dissociative rxns) or unimolecular (units 1/s)
                    if n_reactants == 2:
                        A_old = A_old*1e3
                elif ck_data["rtype"] == "stick":
                    A_old = float(A_old)
                else:
                    raise Exception("reaction type not recognized")

                # units j/mol in chemkin, j/kmol in cantera
                Ea_old = Ea_old*1e3

                # print("new rxn rate", rate)
                a_match = np.isclose(A_old, A_i, rtol=1e-3)
                e_match = np.isclose(Ea_old, Ea_i, rtol=1e-3)

                if a_match and e_match:
                    ck_matches[num] = ck_num
                else: 
                    ck_no_match.append(ck_items)                    

            elif rev: 
                raise Exception("reverse of chemkin reaction found in ct")
            else:
                raise Exception("chemkin reaction does not match ct")

        no_match_ck = list(set(range(len(self.sbr.surf.reactions()))) - set(ck_matches.keys()))
        no_match_ct = list(set(range(len(ck_rule_dict))) - set(ck_matches.values()))
        no_match_dict = dict(zip(no_match_ct, no_match_ck))

        for key, value in no_match_dict.items():
            rxstr = list(ck_rule_dict[value].keys())[0]
            rxstr_ct = self.sbr.surf.reactions()[key].equation
            n_reactants = len(ck_rule_dict[value][rxstr]["reactants"])
            ck_A = ck_rule_dict[value][rxstr]["reac_dict"]["A"]
            ck_Ea = ck_rule_dict[value][rxstr]["reac_dict"]["Ea"]
            ck_dHrxn = ck_rule_dict[value][rxstr]["dHrxn"]
            ct_A = self.sbr.surf.reactions()[key].rate.pre_exponential_factor
            ct_Ea = self.sbr.surf.reactions()[key].rate.activation_energy
            ct_dHrxn = self.sbr.surf.delta_standard_enthalpy[key]

            ck_Ea = ck_Ea*1e3
            if ct_A > 1.0: 
                ck_A = 10**ck_A
                if n_reactants ==2:
                    ck_A = ck_A*1e3

            ck_dHrxn = ck_dHrxn*1e3
            rxrstr_ck = rxstr.replace(" ", "").replace("<=>", "=")
            show_enth = True
            A_tol = 1e-3
            Ea_tol = 5e-2 # 
            if not np.isclose(ck_A, ct_A, rtol=A_tol) and not np.isclose(ck_Ea, ct_Ea, rtol=Ea_tol):
                print(f"rxn {key}, {value}: {rxrstr_ck }, {rxstr_ct}")
                print("both A and Ea don't match")
                print(f"ck: {ck_A}, {ck_Ea}")
                print(f"ct: {ct_A}, {ct_Ea}/n")
                if not np.isclose(ck_dHrxn, ct_dHrxn, rtol=1e-3) or show_enth:
                    print("dHrxn doesn't match")
                    print(f"ck dHrxn: {ck_dHrxn}")
                    print(f"ct dHrxn: {ct_dHrxn}\n")
                    
            elif not np.isclose(ck_A, ct_A, rtol=A_tol) and np.isclose(ck_Ea, ct_Ea, rtol=Ea_tol):
                print(f"rxn {key} {value}: {rxrstr_ck}, {rxstr_ct}")
                print("A doesn't match")
                print(f"ck: {ck_A}")
                print(f"ct: {ct_A}/n")

            elif np.isclose(ck_A, ct_A, rtol=A_tol) and not np.isclose(ck_Ea, ct_Ea, rtol=Ea_tol):
                print(f"rxn {key} {value}: {rxrstr_ck}, {rxstr_ct}")
                print("Ea doesn't match")
                print(f"ck: {ck_Ea/1e6} kj/mol")
                print(f"ct: {ct_Ea/1e6} kj/mol")
                print(f"diff: {abs((ck_Ea-ct_Ea)/1e6)} kj/mol")
                if not np.isclose(ck_dHrxn, ct_dHrxn, rtol=1e-3) or show_enth:
                    print("dHrxn values")
                    print(f"ck dHrxn: {ck_dHrxn/1e6} kj/mol")
                    print(f"ct dHrxn: {ct_dHrxn/1e6} kj/mol\n")

        print(len(ck_rule_dict), len(self.sbr.surf.reactions()))
        # assert len(ck_matches) == len(self.sbr.surf.reactions()), f"not all reactions matched , length of matches: {len(ck_matches)} , length of reactions: {len(self.sbr.surf.reactions())}"
        return ck_matches, ck_no_match