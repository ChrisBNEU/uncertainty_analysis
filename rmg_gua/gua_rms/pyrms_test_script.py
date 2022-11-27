#!/usr/bin/env python
# coding: utf-8

# ## PyRMS test script

# simple test script for executing rms from python, for use in the uncertainty pipeline. rms does sensitivities faster

from pyrms import rms
from diffeqpy import de
from julia import Main
import yaml
from julia import Sundials
from diffeqpy import de
import time 
import matplotlib


class sbr: 
    """
    A minimal stirred batch reactor
    Takes away some of the options from the regular sbr class
    """
    def __init__(
        self,
        rms_file,
        conditions, 
        rtol=1.0e-11,
        atol=1.0e-22,
        time=600,
        ):
        
        # initialize phase dictionary 
        # like ct.solution and ct.Interface in cantera
        phase_dict = rms.readinput(rms_file)


        # convert volume flow to molar flow
        conditions["volume_flowrate"]
        FC_temp = 293.15
        conditions["molar_flow"] = conditions["volume_flowrate"] * 1.01325e5 / (8.3145 * FC_temp) 
        
        # for pyrms need the interface key
        intfc_key = list(phase_dict.keys())[0]

        # mechanism dictionaries index:  
        # phase_dict[phasename]["Species" or "Reactions"]
        gasspcs = phase_dict["gas"]["Species"] 
        gasrxns = phase_dict["gas"]["Reactions"]
        surfacespcs = phase_dict["surface"]["Species"]
        surfacerxns = phase_dict["surface"]["Reactions"]
        interfacerxns = phase_dict[intfc_key]["Reactions"]

        # Define the phase (how species thermodynamic and kinetic properties calculated)
        ig = rms.IdealGas(gasspcs,gasrxns,name="gas") 
        cat = rms.IdealSurface(surfacespcs, surfacerxns, 2.943e-5, name="surface")

        # Set simulation Initial Temp and Pressure
        initialcondsgas = {
                "T":conditions["temperature"],
                "P":conditions["pressure"],
                "CO":conditions["species"]["CO"],
                "CO2":conditions["species"]["CO2"],
                "H2":conditions["species"]["H2"],
        }
        # Define the domain (encodes how system thermodynamic properties calculated)
        domaingas,y0gas,pgas = rms.ConstantTPDomain(phase=ig,initialconds=initialcondsgas,sensitivity=True)

        
        V = conditions["volume"]
        A = conditions["catalyst_area"]
        initialconds = {
                "T":conditions["temperature"],
                "A":conditions["catalyst_area"],
                "X":cat.sitedensity*A
        } #Set simulation Initial Temp and Pressure
        # Define the domain (encodes how system thermodynamic properties calculated)
        domaincat,y0cat,pcat = rms.ConstantTAPhiDomain(phase=cat,initialconds=initialconds,sensitivity=True);


        # ## make reactor, inlet and outlet
        # - makes an anonymous function x->42, is that velocity in? need to check if it is velocity or volume flowrate
        # - also, I think the ```phi``` refers to chemical potential, but I should check, I think constantTPhi is just const T for our case. 

        # In[10]:


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
        t1 = time.time()
        sol = de.solve(react.ode,de.CVODE_BDF(),abstol=1e-20,reltol=1e-8)
        t2 = time.time()


        # In[13]:


        ssys = rms.SystemSimulation(sol,domains,interfaces,p)


        # In[14]:


        rms.plotmolefractions(ssys.sims[1],tol=0.001)


        # In[15]:


        rms.plotrxntransitorysensitivities(ssys, "CC", 100)


        # we need to extract the reaction objects from julia. the sensitivities are just the second index in the tuple, but the first index is the julia reaction objects. I need to spit out a reaction equation for the spreadsheet. during postprocessing, I need to get the sensitivities and trace them back to their originating family

        # In[16]:


        senss = rms.getrxntransitorysensitivities(ssys, "CC", 100)


        # In[26]:


        len(senss)


        # In[28]:


        for index, rxn in enumerate(senss[0]): 
            print(rms.getrxnstr(rxn))


        # In[ ]:




        ## test sbr
        file_dir = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/rms/chem53.rms"