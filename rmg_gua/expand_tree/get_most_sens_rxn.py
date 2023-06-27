#!/usr/bin/env python
# coding: utf-8

# ## PyRMS test script

# simple test script for executing rms from python, for use in the uncertainty pipeline. rms does sensitivities faster
# 
from pyrms import rms
from diffeqpy import de
from julia import Main
import yaml
from julia import Sundials
from diffeqpy import de
import time 
import glob
import matplotlib
from copy import deepcopy
from rmgpy.chemkin import load_chemkin_file
from rmgpy.rmg.model import ReactionModel
from rmgpy.species import Species
from rmgpy.kinetics import StickingCoefficientBEP, StickingCoefficient, SurfaceArrheniusBEP, SurfaceArrhenius
from rmgpy.data.kinetics.database import KineticsDatabase
import os
import copy
import pickle

import sys
repo_dir = os.path.dirname(os.path.dirname(os.path.abspath("")))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(""))))
from rmg_gua.gua_peuqse.peuqse_utilities import get_highest_rms_file

def get_rxn_sensitivities(rms_path, rxn_index):
    if not rms_path: 
        rms_path = os.path.join(repo_dir, "rmg_gua", "baseline")
    else:
        rms_path = os.path.join(repo_dir, "rmg_gua", rms_path)

    file_dir = get_highest_rms_file(rms_path)
    phase_dict = rms.readinput(file_dir)

    base_path = "/".join(file_dir.split("/")[:-2])
    base_path = os.path.join(base_path, "chemkin")
    cmkn_path = os.path.join(base_path, "chem_annotated-gas.inp")
    cmkn_surf_path = os.path.join(base_path, "chem_annotated-surface.inp")
    cmkn_dict_path = os.path.join(base_path, "species_dictionary.txt")
    # file_dir.replace(file_dir.split("/")[-1], "/chemkin")

    expt_condts = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/gua_cantera/all_experiments_reorg_sbr.yaml"

    with open(expt_condts, 'r') as file:
        data = yaml.safe_load(file)

    # pick just one experiment for example 
    conditions = data[2]

    # convert volume flow to molar flow
    conditions["volume_flowrate"]

    FC_temp = 293.15
    conditions["molar_flow"] = conditions["volume_flowrate"] * 1.01325e5 / (8.3145 * FC_temp) 

    intfc_key = list(phase_dict.keys())[0]

    # mechanism dictionaries index:  phase_dict[phasename]["Species" or "Reactions"]
    gasspcs = phase_dict["gas"]["Species"]
    gasrxns = phase_dict["gas"]["Reactions"]
    surfacespcs = phase_dict["surface"]["Species"]
    surfacerxns = phase_dict["surface"]["Reactions"]
    interfacerxns = phase_dict[intfc_key]["Reactions"]

    #Define the phase (how species thermodynamic and kinetic properties calculated)
    ig = rms.IdealGas(gasspcs,gasrxns,name="gas") 
    cat = rms.IdealSurface(surfacespcs, surfacerxns, 2.943e-5, name="surface")

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


    # ## make reactor, inlet and outlet
    # - makes an anonymous function x->42, is that velocity in? need to check if it is velocity or volume flowrate
    # - also, I think the ```phi``` refers to chemical potential, but I should check, I think constantTPhi is just const T for our case. 
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

    # Define domains and interfaces
    domains = (domaingas,domaincat)
    interfaces = [inter,inletgas,outletgas]

    # create a reactor for the system
    react,y0,p = rms.Reactor(domains,(y0gas,y0cat),(0.0,100),interfaces,(pgas,pcat,pinter)) # Create the reactor object

    # run the simulation
    t1 = time.time()
    sol = de.solve(react.ode,de.CVODE_BDF(),abstol=1e-20,reltol=1e-8)
    t2 = time.time()
    print("elapsed time for sim: ", t2-t1)

    ssys = rms.SystemSimulation(sol,domains,interfaces,p)

    rms.getfluxdiagram(ssys,100)


    # load the chemkin file for the mechanism
    chemkin_file = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/chemkin/chem_annotated-gas.inp"
    chemkin_surf_file = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/chemkin/chem_annotated-surface.inp"
    chemkin_dict = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/chemkin/species_dictionary.txt"
    model = ReactionModel()

    model.species, model.reactions = load_chemkin_file(
        chemkin_file,
        chemkin_dict, 
        surface_path=chemkin_surf_file,
        )

    def make_spc(spc): 
        """
        make an RMG object from the rms object
        
        """
        if len(spc.adjlist) > 0:
            rmg_spc = Species().from_adjacency_list(spc.adjlist)
        else: 
            rmg_spc = Species().from_smiles(spc.smiles)

        return rmg_spc

    sens_times = [1e-2, 3e-2, 100]
    sens_rxn_dict = {}
    max_len = 0
    # get the ordering of reactions that are sensitive to CH3OH at different time scales. 
    # our ordering is from highest sensitivity to low. so, start at one. if it is 
    # found to be the first most sensitive at the first time value, then the 5th most sensitive at the 
    # second time value, then the tenth at the third time value then the score it gets is 
    # 1*5*10 = 50

    for time in sens_times:
        sens_rxns, rxn_sens = rms.getrxntransitorysensitivities(ssys, "CH3OH", time, tol = 0)
        counter = 1
        for (rxn, sens) in zip(sens_rxns, rxn_sens):
            if rms.getrxnstr(rxn) in sens_rxn_dict.keys():
                old_sens = sens_rxn_dict[rms.getrxnstr(rxn)][0]
                sens_rxn_dict[rms.getrxnstr(rxn)][0] = counter*old_sens
                sens_rxn_dict[rms.getrxnstr(rxn)][1].append(sens)
                sens_rxn_dict[rms.getrxnstr(rxn)][2].append(counter)
            else:
                # print("reac str: ", rms.getrxnstr(rxn))
                # rxn_adj = rms.getrxnadjlist(rxn)
                # rxn_smiles = rms.getrxnsmiles(rxn)
                # print(rxn_smiles)
                reac_spec = [make_spc(reac) for reac in rxn.reactants]
                prod_spec = [make_spc(prod) for prod in rxn.products]
                sens_rxn_dict[rms.getrxnstr(rxn)] = [counter, [sens],[counter], reac_spec, prod_spec]

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

    min_key = min(sens_cmkn_dict, key=sens_cmkn_dict.get)
    most_sens_rxn = sens_cmkn_dict[min_key][1]
    
    return most_sens_rxn

