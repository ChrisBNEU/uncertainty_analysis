# minimal stirred batch reactor to analyze 5000 rmg runs
###############################################
# Stirred batch reactor script
# Chris Blais, Sevy Harris
# Northeastern University
#
# runs a spinning basket reactor to the specifications in:
# Kinetics of low-pressure methanol synthesis
# Gh Graaf; Ej Stamhuis; Acm Beenackers
# 1988
# 10.1016/0009-2509(88)85127-3
###############################################


import cantera as ct
import math
import yaml
import numpy as np
import os
import sys
from copy import deepcopy
import collections
class MinSBR:
    """
    A minimal stirred batch reactor
    Takes away some of the options from the regular sbr class
    """
    def __init__(
        self,
        yaml_file,
        reac_config, 
        rtol=1.0e-11,
        atol=1.0e-22,
        reactor_type=1,
        energy="off",
        time=600,
        reaction_list = None, 
        results_path = None, 
    ):
        """
        initialize sbr object
        yaml_file = cti or yaml file for mechanism
        reac_config = yaml file containing experimental configuration values
                      Temp, Pressure, concentrations, etc. 
        new_rate_dict = dict of new rates that we want to set
        """

        # load the cantera mechanism
        prefix = os.path.dirname(os.path.abspath(__file__))
        if os.path.exists("/work"):
            self.prefix = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
        else:
            self.prefix = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/"

        if results_path:
            self.results_path = results_path
        else:
            self.results_path = os.getcwd()
        self.expt_id = reac_config['expt_name']
        self.temperature = reac_config['temperature']
        self.pressure = reac_config['pressure']
        self.volume_flow = reac_config['volume_flowrate']
        self.use_for_opt = reac_config['use_for_opt']

        # can probably generalize with isomorphism
        # do that later
        self.x_H2 = reac_config['species']['H2']
        self.x_CO2 = reac_config['species']['CO2']
        self.x_CO = reac_config['species']['CO']
        if 'H2O' in reac_config['species'].keys():
            self.x_H2O = reac_config['species']['H2O']
        else: 
            self.x_H2O = 0.0
        
        # load the experimental TOFs
        self.graaf_meoh_tof = reac_config['output']['CH3OH']
        self.graaf_h2o_tof = reac_config['output']['H2O']

        # molecular weights for mass flow calculations
        MW_CO = 28.01e-3  # [kg/mol]
        MW_CO2 = 44.01e-3  # [kg/mol]
        MW_H2 = 2.016e-3  # [kg/mol]
        MW_H2O = 18.01528e-3  # [kg/mol]

        # define ratios of reactants for plots
        # CO2/(CO+CO2) and (CO2+CO)/H2
        self.CO2_ratio = self.x_CO2 / (self.x_CO + self.x_CO2)
        self.H2_ratio = (self.x_CO2 + self.x_CO) / self.x_H2

        # create thermo phases
        self.yaml_file = yaml_file
        self.gas = ct.Solution(yaml_file, "gas")
        self.surf = ct.Interface(yaml_file, "surface1", [self.gas])
        
        # modify the reactions if specified
        if reaction_list is not None:
            self.gas.TP = 298.0, 101325.0
            self.surf.TP = 298.0, 101325.0
            rxn_pert_dict, thermo_pert_dict = self.load_perts(reaction_list)
            self.change_be(thermo_pert_dict)
            self.change_reactions_rmg(rule_dict=rxn_pert_dict)

        # pull out species names
        for spec_str in self.gas.species_names:
            if spec_str.startswith("CO("):
                self.co_str = spec_str
            if spec_str.startswith("CO2("):
                self.co2_str = spec_str
            if spec_str.startswith("H2("):
                self.h2_str = spec_str
            if spec_str.startswith("H2O("):
                self.h2o_str = spec_str
            if spec_str.startswith("CH3OH("):
                self.ch3oh_str = spec_str
        
        # CO/CO2/H2/H2: typical is
        self.concentrations_rmg = {
            self.co_str: self.x_CO,
            self.co2_str: self.x_CO2,
            self.h2_str: self.x_H2,
            self.h2o_str: self.x_H2O,
        }
        # initialize T and P
        self.gas.TPX = self.temperature, self.pressure, self.concentrations_rmg
        self.surf.TP = self.temperature, self.pressure

        # if a mistake is made with the input moles,
        # cantera will normalize the mole fractions.
        # make sure that we are reporting
        # the normalized values
        self.x_CO = float(self.gas[self.co_str].X)
        self.x_CO2 = float(self.gas[self.co2_str].X)
        self.x_H2 = float(self.gas[self.h2_str].X)
        self.x_H2O = float(self.gas[self.h2o_str].X)

        # create gas inlet
        self.inlet = ct.Reservoir(self.gas)

        # create gas outlet
        self.exhaust = ct.Reservoir(self.gas)

        # Reactor volume
        self.rvol = reac_config['volume'] 

        # Catalyst Surface Area
        self.cat_area = reac_config['catalyst_area']  # [m^3]
        self.cat_area_str = "%s" % "%.3g" % self.cat_area
        self.total_sites = (self.cat_area*self.surf.site_density)*1e3 #[mol sites]

        # reactor initialization
        # always use an IdealGasReactor
        # TODO check about reactor type, see if you can get rid of this block
        self.reactor_type = reactor_type
        self.energy = energy
        if self.reactor_type == 0:
            self.r = ct.Reactor(self.gas, energy=self.energy)
            self.reactor_type_str = "Reactor"
        elif reactor_type == 1:
            self.r = ct.IdealGasReactor(self.gas, energy=self.energy)
            self.reactor_type_str = "IdealGasReactor"
        elif reactor_type == 2:
            self.r = ct.ConstPressureReactor(self.gas, energy=self.energy)
            self.reactor_type_str = "ConstPressureReactor"
        elif reactor_type == 3:
            self.r = ct.IdealGasConstPressureReactor(self.gas, energy=self.energy)
            self.reactor_type_str = "IdealGasConstPressureReactor"

        # calculate the available catalyst area in a differential reactor
        self.rsurf = ct.ReactorSurface(self.surf, self.r, A=self.cat_area)
        self.r.volume = self.rvol

        self.surf.coverages = "X(1):1.0"

        # flow controllers (Graaf measured flow at 293.15 and 1 atm)
        FC_temp = 293.15
        self.molar_flow = self.volume_flow * ct.one_atm / (8.3145 * FC_temp)  # [mol/s]
        self.mass_flow = self.molar_flow * (
            self.x_CO * MW_CO + self.x_CO2 * MW_CO2 + self.x_H2 * MW_H2 + self.x_H2O * MW_H2O
        )  # [kg/s]
        self.mfc = ct.MassFlowController(self.inlet, self.r, mdot=self.mass_flow)

        # A PressureController has a baseline mass flow rate matching the 'master'
        # MassFlowController, with an additional pressure-dependent term. By explicitly
        # including the upstream mass flow rate, the pressure is kept constant without
        # needing to use a large value for 'K', which can introduce undesired stiffness.
        self.outlet_mfc = ct.PressureController(self.r, self.exhaust, master=self.mfc, K=0.01)

        # initialize reactor network
        self.sim = ct.ReactorNet([self.r])

        # set relative and absolute tolerances on the simulation
        self.sim.rtol = rtol
        self.sim.atol = atol
        
        self.time=time


    def load_perts(self, pert_list):
        """ 
        load the perturbations to a dictionary
        """
        ck_rule_dict_path = os.path.join(self.results_path, "ck_rule_dict.yaml")
        with open(ck_rule_dict_path, 'r') as f:
            ck_rule_dict = yaml.load(f, Loader=yaml.FullLoader)

        rule_config_file = os.path.join(self.results_path, "rule_config.yaml")
        with open(rule_config_file, "r") as f:
            rule_dict_orig = yaml.safe_load(f)

        # # load in the perturbations. haven't thought of a more intelligent way to do this, 
        # but order is preserved in yaml for dict for python 3.6+ if we safe load 
        # and save with sort_keys = False
        rule_dict = deepcopy(rule_dict_orig)
        ck_matches = {}
        ck_no_match = []
        count = 0

        if len(pert_list) > 0:
            assert len(pert_list) == len(rule_dict_orig)*3 + 5, "number of reactions in rxn_list does not match number of reactions in rule_dict_orig"
            print("changing reactions")
            for rule, vals in rule_dict_orig.items():
                rule_dict[rule]["A"] = pert_list[count]
                count += 1
                rule_dict[rule]["E0"] = pert_list[count]
                count += 1
                rule_dict[rule]["alpha"] = pert_list[count]
                count += 1
            # load the binding energy perturbations last values in dict
            thermo_pert_dict = {'C': 0., 'O': 0., 'N': 0., 'H': 0., 'vdw': 0.}
            for num, (key, value) in enumerate(thermo_pert_dict.items()):
                thermo_pert_dict[key] = pert_list[count]
                count += 1


        return rule_dict, thermo_pert_dict

    def change_be(self, thermo_pert_dict):
        """
        change the species BE by altering the enthalpy of formation
        """
        # load in the test species yaml file if requested
        test_spec_path = os.path.join(self.results_path, "be_test_config.yaml")        
        with open(test_spec_path, "r") as f:
            test_spec = yaml.safe_load(f)

        # load the lookup for which species have which bond order
        thermo_bo_dict_path = os.path.join(self.results_path, "thermo_pert.yaml")
        with open(thermo_bo_dict_path, 'r') as f:
            thermo_bo_dict = yaml.load(f, Loader=yaml.FullLoader)
        
        # dh is range of 0.3 eV, or 1e7 j/kmol
        kj_thermo_pert = {k: v/1e6 for k, v in thermo_pert_dict.items()}

        old_h298 = {}
        for spec in self.surf.species():
            if spec.name == "X(1)":
                print(spec.name)
                pass
            else:
                dh = 0. 
                species = spec.name
                if test_spec and species in test_spec.values():
                    old_h298[species] = spec.thermo.h(298.15)/1e6
                # calculate dh
                bo_spec = thermo_bo_dict[species] 
                for atom, value in bo_spec.items():
                    dh += thermo_pert_dict[atom] * value

                st_orig = spec.thermo

                coeffs = st_orig.coeffs
                coeffs[[6, 13]] += dh / ct.gas_constant
                s_new = ct.NasaPoly2(st_orig.min_temp, st_orig.max_temp, st_orig.reference_pressure, coeffs)
                spec.thermo = s_new
                self.surf.modify_species(self.surf.species_index(species), spec)

                if test_spec and species in test_spec.values():
                    diff = spec.thermo.h(298.15)/1e6 - old_h298[species]
                    pert = dh/1e6
                    assert np.isclose(diff, pert, rtol=1e-3, atol=1e-3)
        
        return old_h298

    def scale_ea(self, E0, alpha, rxn, rxind):
        """ 
        scale the activation energy for each reaction using either blowers masel
        or BEP
        """
        hrxn = self.surf.delta_standard_enthalpy[rxind]
        reactants = rxn.reactants
        products = rxn.products
        # get enthalpy of reaction at 0K
        reactants = rxn.reactants
        products = rxn.products
        reac_enth = 0
        for reactant, stoic in reactants.items():
            try:
                ind = self.surf.species_index(reactant)
                reac_enth += self.surf.species(ind).thermo.h(0.01) * stoic
            except ValueError:
                ind = self.gas.species_index(reactant)
                reac_enth += self.gas.species(ind).thermo.h(0.01) * stoic
        prod_enth = 0
        for product, stoic in products.items():
            try:
                ind = self.surf.species_index(product)
                prod_enth += self.surf.species(ind).thermo.h(0.01) * stoic
            except ValueError:
                ind = self.gas.species_index(product)
                prod_enth += self.gas.species(ind).thermo.h(0.01) * stoic

        H0 = prod_enth - reac_enth

        # scale Ea
        Ea = E0 + alpha * hrxn

        # check if Ea is negative or less than rxn enthalpy at 0K
        if H0 > 0.0 and Ea < H0:
            Ea = H0
        elif H0 < 0.0 and Ea < 0.0:
            Ea = 0
        return Ea

    def change_reactions_rmg(self, rule_dict = [], test = False):
        """
        find the template for a given reaction
        for now, we don't really care about gas phase rxns
        """
        # load the chenkin rxn_config file
        ck_rule_dict_path = os.path.join(self.results_path, "ck_rule_dict.yaml")
        with open(ck_rule_dict_path, 'r') as f:
            ck_rule_dict = yaml.load(f, Loader=yaml.FullLoader)

        ck_matches = {}
        ck_no_match = []
        for ck_num, ck_items in enumerate(ck_rule_dict):
            ck_rxn = list(ck_items.keys())[0] 
            ck_data = list(ck_items.values())[0]
            ck_reac = collections.Counter(ck_data["reactants"])
            ck_prod = collections.Counter(ck_data["products"])

            # for num, rxn in enumerate(self.surf.reactions()):
            num = ck_num
            rxn = self.surf.reactions()[ck_num]
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

            # if we have a match, write the A, Ea to the ct rxn
            if fwd:
                if num in ck_matches.keys():
                    print(f"duplicate reaction found {ck_rxn}, moving to next rxn")
                    pass
                else: 
                    new_rxn = self.surf.reactions()[num]
                    # print("oldrxn: ", self.surf.reactions()[num].rate)

                    A_old = self.surf.reactions()[num].rate.pre_exponential_factor
                    Ea_old = self.surf.reactions()[num].rate.activation_energy

                    # pull the A, Ea from the chemkin file
                    A_src = rule_dict[ck_data["source"]]["A"]
                    E0_src = rule_dict[ck_data["source"]]["E0"]
                    alpha_src = rule_dict[ck_data["source"]]["alpha"]
                    # hrxn = ck_data["dHrxn"]*1e3
                    hrxn = self.surf.delta_standard_enthalpy[num] 

                    # get nreactants for A unit consistency
                    n_reactants = len(ck_data["reactants"])

                    # record initial values from cantera
                    A_i = self.surf.reactions()[num].rate.pre_exponential_factor
                    Ea_i = self.surf.reactions()[num].rate.activation_energy
                    b_i = self.surf.reactions()[num].rate.temperature_exponent
                    
                    if ck_data["rtype"] == "arrhenius":
                        A_i = 10**float(A_src)
                        # for bimolecular reactions we have m^2/mol*s in chemkin, 
                        # but m^2/kmol*s in cantera. all others are either sticking 
                        # (trimolecular dissociative rxns) or unimolecular (units 1/s)
                        if n_reactants == 2:
                            A_i = A_i*1e3
                    elif ck_data["rtype"] == "stick":
                        A_i = float(A_src)
                    else:
                        raise Exception("reaction type not recognized")
                    
                    # account for rxn path degeneracy
                    if ck_data["deg"] > 1.0:
                        A_i = A_i * ck_data["deg"]
                    
                    # units j/mol in chemkin, j/kmol in cantera
                    E0_src = E0_src*1e3
                    Ea_i = self.scale_ea(E0_src, alpha_src, rxn, num)

                    rate = ct.Arrhenius(A = A_i, E = Ea_i, b = b_i)
                    new_rxn.rate = rate
                    self.surf.modify_reaction(num, new_rxn)

                    # print("new rxn rate", rate)
                    a_match = np.isclose(A_old, A_i, rtol=1e-3)
                    e_match = np.isclose(Ea_old, Ea_i, rtol=1e-3)

                    if a_match and e_match:
                        ck_matches[num] = ck_num
                    else: 
                        ck_no_match.append(ck_items)   
                        ck_matches[num] = ck_num

            elif rev: 
                raise Exception("reverse of chemkin reaction found in ct")
            else:
                raise Exception("chemkin reaction does not match ct")

        return ck_matches


    def run_reactor_ss_memory(self):
        """
        Run single reactor to steady state and save the results to an ordered dictionary in memory
        """

        # how to sort out gas and surface such that we can attach units?
        gas_ROP_str = [i + " ROP [kmol/m^3 s]" for i in self.gas.species_names]

        # Okay, this is weird. rsurf.kinetics.net_production_rates includes both gas and 
        # surface rates, but surf.species_names only has surface species
        all_surf_rop_specs = self.gas.species_names + self.surf.species_names
        surf_ROP_str = [i + " ROP [kmol/m^2 s]" for i in all_surf_rop_specs]

        gasrxn_ROP_str = [i + " ROP [kmol/m^3 s]" for i in self.gas.reaction_equations()]
        surfrxn_ROP_str = [i + " ROP [kmol/m^2 s]" for i in self.surf.reaction_equations()]

        # run the simulation
        # 600s chosen because residence time is ~20 seconds
        # reactor volume is 1.35*10^-4 m^3
        # volume flowrate is always greater than 1e-6 m^3/s, meaning at most 
        # our residence time is 134s. so 134*4 = 536.  
        self.sim.advance(self.time)
        results = {}
        results['experiment'] = self.expt_id
        results['use_for_opt'] = self.use_for_opt
        results['time (s)'] = self.sim.time
        results['T (K)'] = self.temperature
        results['P (Pa)'] = self.gas.P
        results['V (m^3/s)'] = self.volume_flow
        results['x_CO initial'] = self.x_CO
        results['x_CO2 initial'] = self.x_CO2
        results['x_H2 initial'] = self.x_H2
        results['x_H2O initial'] = self.x_H2O
        results['CO2/(CO2+CO)'] = self.CO2_ratio
        results['(CO+CO2/H2)'] = self.H2_ratio
        results['T (K) final'] = self.gas.T
        results['Rtol'] = self.sim.rtol
        results['Atol'] = self.sim.atol
        results['reactor type'] = self.reactor_type_str
        results['energy on?'] = self.energy
        results['catalyst area'] = self.cat_area
        results['graaf MeOH TOF 1/s'] = self.graaf_meoh_tof 
        results['graaf H2O TOF 1/s'] = self.graaf_h2o_tof
        
        results['RMG MeOH TOF 1/s'] = float(self.molar_flow*(self.r.thermo[self.ch3oh_str].X)/(self.total_sites))
        results['RMG H2O TOF 1/s'] = float(self.molar_flow*(self.r.thermo[self.h2o_str].X - self.x_H2O)/(self.total_sites))
        
        # results['RMG MeOH TOF 1/s'] = self.rsurf.kinetics.net_production_rates[self.gas.species_index(self.ch3oh_str)]/self.surf.site_density
        # results['RMG H2O TOF 1/s'] = self.rsurf.kinetics.net_production_rates[self.gas.species_index(self.h2o_str)]/self.surf.site_density
        
        
        results['error squared MeOH TOF'] = ((results['graaf MeOH TOF 1/s'] - results['RMG MeOH TOF 1/s'])/results['graaf MeOH TOF 1/s'] )**2
        results['error squared H2O TOF'] = ((results['graaf H2O TOF 1/s'] - results['RMG H2O TOF 1/s'])/results['graaf H2O TOF 1/s'])**2
        results['obj_func'] = results['error squared MeOH TOF'] + results['error squared H2O TOF']

        # adding alternative objective function using avg of log(rate_rmg/rate_exp)
        results['log10(RMG/graaf) MeOH TOF'] = np.log10(max(1e-9, results['RMG MeOH TOF 1/s']/results['graaf MeOH TOF 1/s']))
        results['log10(RMG/graaf) H2O TOF'] = np.log10(max(1e-9, results['RMG H2O TOF 1/s']/results['graaf H2O TOF 1/s']))
        results['log10(RMG/graaf) TOF'] = 0.5 * ( results['log10(RMG/graaf) MeOH TOF'] + results['log10(RMG/graaf) H2O TOF'])
        
        for i in range(0, len(self.gas.X)):
            results[self.gas.species_names[i]] = self.gas.X[i]
        for i in range(0, len(self.surf.X)):
            results[self.surf.species_names[i]] = self.surf.X[i]
        
        # Enter the ROP's
        for i in range(0, len(self.r.kinetics.net_production_rates)):
            results[gas_ROP_str[i]] = self.r.kinetics.net_production_rates[i]

        for i in range(0, len(self.rsurf.kinetics.net_production_rates)):
            results[surf_ROP_str[i]] = self.rsurf.kinetics.net_production_rates[i]

        for i in range(0, len(self.rsurf.kinetics.net_rates_of_progress)):
            results[surfrxn_ROP_str[i]] = self.rsurf.kinetics.net_rates_of_progress[i]

        if gasrxn_ROP_str:
            for i in range(0, len(self.r.kinetics.net_rates_of_progress)):
                results[gasrxn_ROP_str[i]] = self.r.kinetics.net_rates_of_progress[i]

        return results

def run_sbr_test():
    """
    mostly used to run in debugger in vscode
    """
    cti_file_path = "../../../uncertainty_output_folder/run_0001/cantera/chem_annotated.cti"
    
    # load data file and use a random experiment for testing
    with open ('../all_experiments_reorg_sbr.yaml', 'r+') as f:
        y = yaml.safe_load(f)
    

    # initialize reactor
    sbr_ss = MinSBR(
        cti_file_path,
        reac_config = y[1],
        rtol=1.0e-11,
        atol=1.0e-22,
    )

    # run to SS
    results = sbr_ss.run_reactor_ss_memory()
    print("done")

if __name__ == "__main__":
    # execute only if run as a script
    run_sbr_test()
