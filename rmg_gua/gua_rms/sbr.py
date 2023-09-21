from pyrms import rms
from diffeqpy import de
from julia import Main
import yaml
import time 
import matplotlib
import pickle
from copy import deepcopy
import os
import numpy as np
# in peuquse run we cannot use rmgpy objects
from rmgpy.chemkin import load_chemkin_file
from rmgpy.rmg.model import ReactionModel
from rmgpy.species import Species
from rmgpy.kinetics import StickingCoefficientBEP, StickingCoefficient, SurfaceArrheniusBEP, SurfaceArrhenius
from rmgpy.data.kinetics.database import KineticsDatabase

def make_spc(spc): 
    """
    make an RMG object from the rms object

    """
    if len(spc.adjlist) > 0:
        rmg_spc = Species().from_adjacency_list(spc.adjlist)
    else: 
        rmg_spc = Species().from_smiles(spc.smiles)

    return rmg_spc

class rms_sbr:

    def __init__(
        self,
        file_path,
        reac_config,
        rtol=1.0e-11,
        atol=1.0e-22,
    ):

        # establish some values that pass through to output csv: 
        self.file_dir = file_path
        self.expt_id = reac_config['expt_name']
        self.temperature = reac_config['temperature']
        self.pressure = reac_config['pressure']
        self.volume_flow = reac_config['volume_flowrate']
        self.use_for_opt = reac_config['use_for_opt']

        self.reactime = 600 # seconds
        self.reactor_type_str = "IdealGasReactor"
        self.energy = "off"
        # can probably generalize with isomorphism
        # do that later
        self.x_H2 = reac_config['species']['H2']
        self.x_CO2 = reac_config['species']['CO2']
        self.x_CO = reac_config['species']['CO']
        if 'H2O' in reac_config['species'].keys():
            self.x_H2O = reac_config['species']['H2O']
        else: 
            self.x_H2O = 0.0
        
        self.CO2_ratio = self.x_CO2 / (self.x_CO + self.x_CO2)
        self.H2_ratio = (self.x_CO2 + self.x_CO) / self.x_H2
        
        # load the experimental TOFs
        self.graaf_meoh_tof = reac_config['output']['CH3OH']
        self.graaf_h2o_tof = reac_config['output']['H2O']

        # load the experimental mole fracs out
        self.graaf_mole_fracs = reac_config['species_out']

        # make phase dict
        phase_dict = rms.readinput(self.file_dir)

        # convert volume flow to molar flow
        FC_temp = 293.15
        reac_config["molar_flow"] = reac_config["volume_flowrate"] * 1.01325e5 / (8.3145 * FC_temp) 
        self.molar_flow = reac_config["molar_flow"]
        # get rxns and species from gas surf and interface phases
        # mechanism dictionaries index:  
        # phase_dict[phasename]["Species" or "Reactions"]

        # This is a workaround for pyrms, have to feed in interface object key
        intfc_key = list(phase_dict.keys())[0]
        gasspcs = phase_dict["gas"]["Species"]
        gasrxns = phase_dict["gas"]["Reactions"]
        surfacespcs = phase_dict["surface"]["Species"]
        surfacerxns = phase_dict["surface"]["Reactions"]
        interfacerxns = phase_dict[intfc_key]["Reactions"]

        # get site density from chemkin surf file
        self.base_path = "/".join(self.file_dir.split("/")[:-2])
        chemkin_surf_file = os.path.join(
            self.base_path, "chemkin", "chem_annotated-surface.inp")
        
        # open chemkin surf file. first line has site density
        with open(chemkin_surf_file, 'r') as f:
            chemkin_surf_file_lines = f.readlines()
            for line in chemkin_surf_file_lines:
                if line.startswith("SITE"):
                    self.site_density = float(line.split("/")[1])
                    self.site_density = self.site_density * 1e4  # convert from mol/cm2 to mol/m2
                    break
        
        # Define the phase (how species thermodynamic and kinetic properties calculated)
        ig = rms.IdealGas(gasspcs,gasrxns,name="gas") 
        cat = rms.IdealSurface(surfacespcs, surfacerxns,self.site_density, name="surface")

        # Set simulation gas Initial Temp and Pressure
        initialcondsgas = {
                "T":reac_config["temperature"],
                "P":reac_config["pressure"],
                "CO":reac_config["species"]["CO"],
                "CO2":reac_config["species"]["CO2"],
                "H2":reac_config["species"]["H2"],
        } 

        # Define the gas domain (encodes how system thermodynamic properties calculated)
        domaingas,y0gas,pgas = rms.ConstantTPDomain(phase=ig,initialconds=initialcondsgas,sensitivity=True)

        # Set simulation surf Initial Temp and Pressure
        V = reac_config["volume"]
        A = reac_config["catalyst_area"]
        self.cat_area = A
        self.total_sites = (self.cat_area*self.site_density) #[mol sites]
        initialconds = {
                "T":reac_config["temperature"],
                "A":reac_config ["catalyst_area"],
                "X":cat.sitedensity*A
        } 
        # Define the surf domain (encodes how system thermodynamic properties calculated)
        domaincat,y0cat,pcat = rms.ConstantTAPhiDomain(phase=cat,initialconds=initialconds,sensitivity=True);

        
        # ## make reactor, inlet and outlet
        # - makes an anonymous function x->42, is that velocity in? need to check if it is velocity or volume flowrate
        # - also, I think the ```phi``` refers to chemical potential, but I should check, I think constantTPhi is just const T for our case. 

        initialcondsinlet = {
                "T":reac_config["temperature"],
                "P":reac_config["pressure"],
                "CO":reac_config["species"]["CO"],
                "CO2":reac_config["species"]["CO2"],
                "H2":reac_config["species"]["H2"],
            }

        # construct reactor
        inter,pinter = rms.ReactiveInternalInterfaceConstantTPhi(domaingas,domaincat,interfacerxns,initialcondsinlet["T"],A);

        # make inlet and outlet
        inletgas = rms.Inlet(domaingas,initialcondsinlet,Main.eval("x->"+str(reac_config ["molar_flow"])))
        outletgas = rms.Outlet(domaingas,Main.eval("x->"+str(reac_config ["molar_flow"])))

        # Define domains and interfaces
        self.domains = (domaingas,domaincat)
        self.interfaces = [inter,inletgas,outletgas]

        # create a reactor for the system
        self.react,self.y0,self.p = rms.Reactor(
            self.domains,
            (y0gas,y0cat),
            (0.0,self.reactime),
            self.interfaces,
            (pgas,pcat,pinter),
            ) # Create the reactor object
        
        # tolerances
        self.atol = atol
        self.rtol = rtol

    def run_simulation(self, peuqse= False, sens_spcs=["CH3OH", "CC", "CH4"]):
        """
        run the simulation and save the results to a csv, like the cantera script
        peuqse - cuts down on i/o operations so we run quicker
        """
        # run the simulation
        t1 = time.time()
        self.sol = de.solve(
            self.react.ode,
            de.CVODE_BDF(),
            abstol=self.atol,
            reltol=self.rtol,)
        t2 = time.time()
        print("elapsed time for sim: ", t2-t1)

    
        self.ssys = rms.SystemSimulation(self.sol,self.domains,self.interfaces,self.p)
        results = {}
        # get mole fractions from self.ssys for csv
        results['experiment'] = self.expt_id
        results['use_for_opt'] = self.use_for_opt
        results['time (s)'] = self.reactime
        results['T (K)'] = self.temperature
        results['P (Pa)'] = rms.getP(self.ssys.sims[0], self.reactime) # given in Pa
        results['V (m^3/s)'] = self.volume_flow
        results['x_CO initial'] = self.x_CO
        results['x_CO2 initial'] = self.x_CO2
        results['x_H2 initial'] = self.x_H2
        results['x_H2O initial'] = self.x_H2O
        results['CO2/(CO2+CO)'] = self.CO2_ratio
        results['(CO+CO2/H2)'] = self.H2_ratio
        results['T (K) final'] = rms.getT(self.ssys.sims[0], self.reactime)
        results['Rtol'] = self.rtol
        results['Atol'] = self.atol
        results['reactor type'] = self.reactor_type_str
        results['energy on?'] = self.energy
        results['catalyst area'] = self.cat_area
        results['graaf MeOH TOF 1/s'] = self.graaf_meoh_tof 
        results['graaf H2O TOF 1/s'] = self.graaf_h2o_tof    
        meoh_x = rms.molefractions(self.ssys.sims[0],"CH3OH", self.reactime)
        h2o_x = rms.molefractions(self.ssys.sims[0],"H2O", self.reactime)
        results['RMG MeOH TOF 1/s'] = float(self.molar_flow*(meoh_x)/(self.total_sites))
        results['RMG H2O TOF 1/s'] = float(self.molar_flow*(h2o_x- self.x_H2O)/(self.total_sites))

        # objective function: make RMSE of Mole fractions for all species, weighted by their value 
        results['error squared MeOH TOF'] = ((results['graaf MeOH TOF 1/s'] - results['RMG MeOH TOF 1/s'])/results['graaf MeOH TOF 1/s'] )**2
        results['error squared H2O TOF'] = ((results['graaf H2O TOF 1/s'] - results['RMG H2O TOF 1/s'])/results['graaf H2O TOF 1/s'])**2
        results['obj_func'] = results['error squared MeOH TOF'] + results['error squared H2O TOF']

        # adding alternative objective function using avg of log(rate_rmg/rate_exp)
        results['log10(RMG/graaf) MeOH TOF'] = np.log10(max(1e-9, results['RMG MeOH TOF 1/s']/results['graaf MeOH TOF 1/s']))
        results['log10(RMG/graaf) H2O TOF'] = np.log10(max(1e-9, results['RMG H2O TOF 1/s']/results['graaf H2O TOF 1/s']))
        results['log10(RMG/graaf) TOF'] = 0.5 * ( results['log10(RMG/graaf) MeOH TOF'] + results['log10(RMG/graaf) H2O TOF'])

        
        # get mole fractions
        # gas
        gas_names = self.ssys.sims[0].names
        gas_moles = rms.molefractions(self.ssys.sims[0], self.reactime)
        for (name, moles) in zip(gas_names, gas_moles):
            results[name] = moles
        
        # surface
        surf_names = self.ssys.sims[1].names
        surf_moles = rms.molefractions(self.ssys.sims[1], self.reactime)
        for (name, moles) in zip(surf_names, surf_moles):
            results[name] = moles
        
        # get graaf mole fractions
        for spec, val in self.graaf_mole_fracs.items():
            results['graaf moles ' + spec] = val

        species_err = list(self.graaf_mole_fracs.keys())

        # get percent error
        for spec in species_err:
            results["Error % " + spec] = 100*abs(results["graaf moles " + spec] -
                                                 results[spec])/results["graaf moles " + spec]
        results["Sum Error %"] = sum(
            [results["Error % " + spec] for spec in species_err])
         
        # get reaction for CH3OH, ethane, and CH4 at reactor outlet
        rxn_strs = []
        for rxn in self.ssys.reactions: 
            rxn_strs.append(rms.getrxnstr(rxn))

        for spec in sens_spcs:
            # get the sensitivities. transitory sensitivities returns the full sensitivity matrix
            # calculated from the jacobian. returns a (n_species, n_reactions) x (n_species, n_reactions) matrix, 
            # so to just get the reaction sensitivities we single out a species for out row, then
            # slice out the first n_species columns. adapted from plotrxntransitorysensitivities in rms
            try:
                sens_items, _ = rms.transitorysensitivitiesfulltrapezoidal(self.ssys, 600)
            except RuntimeError as e:
                print("error getting sensitivities in julia")
                continue
            
            try:
                ind = self.ssys.names.index(spec)
            except ValueError:
                print(spec, " not found in species list")
                continue
            sens_items_rxn = sens_items[ind,len(self.ssys.names):]

            for rxn, sens in zip(rxn_strs, sens_items_rxn):
                results[rxn + " sens to " + spec] = sens
        
        if peuqse: 
            # get sensitivities
            sens_times = [1e-2, 3e-1, 600]
            sens_rxn_dict = {}
            max_len = 0
            # get the ordering of reactions that are sensitive to CH3OH at different time scales. 
            # our ordering is from highest sensitivity to low. so, start at one. if it is 
            # found to be the first most sensitive at the first time value, then the 5th most sensitive at the 
            # second time value, then the tenth at the third time value then the score it gets is 

            for stime in sens_times:
                ind = self.ssys.names.index("CH3OH")
                sens_items, _ = rms.transitorysensitivitiesfulltrapezoidal(self.ssys, stime)
                sens_items_rxn = sens_items[ind,len(self.ssys.names):]

                counter = 1
                for (rxn, sens) in zip(sens_rxns, rxn_sens):
                    if rxn in sens_rxn_dict.keys():
                        old_sens = sens_rxn_dict[rxn][0]
                        sens_rxn_dict[rxn][0] = counter*old_sens
                        sens_rxn_dict[rxn][1].append(sens)
                        sens_rxn_dict[rxn][2].append(counter)
                    else:
                        rxn_adj = rms.getrxnadjlist(rxn)
                        reac_spec = [make_spc(reac) for reac in rxn.reactants]
                        prod_spec = [make_spc(prod) for prod in rxn.products]
                        sens_rxn_dict[rxn] = [counter, [sens],[counter], reac_spec, prod_spec]

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
            # save the dictionary to a yaml file:
            with open(os.path.join(self.base_path, "sens_rms_dict.pickle"), 'wb') as file:
                pickle.dump(sens_rxn_dict, file)

            # print the dict: 
            for rxn_str, entry in sens_rxn_dict.items():
                results[rxn_str + "CH3OH sens"] = entry[1][len(sens_times)-1]
                # print(rxn_str, entry[0])

            # load the chemkin file for the mechanism
            chemkin_file = os.path.join(self.base_path, "chemkin", "chem_annotated-gas.inp")
            chemkin_surf_file = os.path.join(self.base_path,"chemkin", "chem_annotated-surface.inp")
            chemkin_dict = os.path.join(
                self.base_path, "chemkin", "species_dictionary.txt")

            # build reaction model
            model = ReactionModel()
            model.species, model.reactions = load_chemkin_file(
                chemkin_file,
                chemkin_dict, 
                surface_path=chemkin_surf_file,
                )

            # now match up the rms reaction sensitivities with the actual rxn in chemkin
            # match species function does forward and reverse. if there are multiple matches, it will
            # print a warning.
            match_list = []
            sens_cmkn_dict = {}
            for rxn_str, entry in sens_rxn_dict.items():
                for rxn in model.reactions:
                    counter = 0
                    if rxn.matches_species(entry[3],entry[4]):
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

            # remove entries that are not estimated from families

            # save the dictionary in a pickle file
            cmkn_pickle_path = os.path.join(self.base_path, "sens_cmkn_dict.pickle")
            with open(cmkn_pickle_path, "wb") as handle:
                pickle.dump(sens_cmkn_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return results
    


