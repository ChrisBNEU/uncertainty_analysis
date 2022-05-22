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
    ):
        """
        initialize sbr object
        yaml_file = cti or yaml file for mechanism
        reac_config = yaml file containing experimental configuration values
                      Temp, Pressure, concentrations, etc. 
        """
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
        self.sim.advance_to_steady_state()
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
        results['RMG MeOH TOF 1/s'] = self.rsurf.kinetics.net_production_rates[self.gas.species_index(self.ch3oh_str)]/self.surf.site_density
        results['RMG H2O TOF 1/s'] = self.rsurf.kinetics.net_production_rates[self.gas.species_index(self.h2o_str)]/self.surf.site_density
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
    cti_file_path = "../../uncertainty_output_folder/run_0001/cantera/chem_annotated.cti"
    
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
