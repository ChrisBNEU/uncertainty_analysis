import sys
import os 
if os.path.exists("/work"):
    prefix = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
else:
    prefix = "/Users/blais.ch/Documents/_01_code/05_Project_repos_Github/meOH_repos/uncertainty_analysis/"
sys.path.append(prefix)

from rmg_gua.gua_peuqse.peuqse_utilities import *
import numpy as np 
import pandas as pd
import yaml

# 
sens_spcs = ["CH3OH"]
make_ct_peuq_input(base_model, overwrite=True, sens_specs=sens_spcs)