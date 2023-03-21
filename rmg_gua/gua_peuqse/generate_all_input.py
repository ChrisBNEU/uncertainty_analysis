import collections
import sys
from rmg_gua.gua_rms.sbr import rms_sbr
import os
from pyrms import rms
from julia import Main
import time 
import yaml
import pandas as pd
import math
import cantera as ct
from copy import deepcopy
from rmgpy.chemkin import load_chemkin_file
from rmgpy.rmg.model import ReactionModel
from rmgpy.species import Species
from rmgpy.kinetics import StickingCoefficientBEP, StickingCoefficient, SurfaceArrheniusBEP, SurfaceArrhenius
from rmgpy.data.kinetics.database import KineticsDatabase
from matplotlib import pyplot as plt
from peuqse_utilities import *


# run the script with the location of the base model specied
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("please specify path to base model")
        sys.exit(1)

    base_model = sys.argv[1]

    sens_spcs = ["CH3OH", "CC", "CH4"]
    make_ct_peuq_input(base_model, overwrite=True, sens_specs=sens_spcs)
    make_ct_expt_file(base_model)




