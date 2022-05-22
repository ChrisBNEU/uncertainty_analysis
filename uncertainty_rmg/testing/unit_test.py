import cantera as ct
import numpy as np
from torch.quasirandom import SobolEngine
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.data.kinetics.database import KineticsDatabase
from IPython.display import Image
import seaborn as sns
import pickle

sns.set_palette('Dark2')

# make sobol map
with open('../../sobol_map.pickle', 'rb') as f:
    sobol_map = pickle.load(f)

DELTA_E0_MAX_J_MOL = 30000
DELTA_E0_MAX_J_MOL_VDW = 20000

# Create the pseudo randoms
N = 20
sobol = SobolEngine(dimension=300, scramble=True, seed=100)
x_sobol = sobol.draw(N)

# set perturb values
pert = 4
perturb = f"000{pert}"

E_0_c = float(DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[pert,sobol_map["C_BE"][0]] * DELTA_E0_MAX_J_MOL)/9.6e4
E_0_o = float(DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[pert,sobol_map["O_BE"][0]] * DELTA_E0_MAX_J_MOL)/9.6e4
E_0_h = float(DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[pert,sobol_map["H_BE"][0]] * DELTA_E0_MAX_J_MOL)/9.6e4
E_0_vdw = float(DELTA_E0_MAX_J_MOL_VDW  - 2.0 * x_sobol[pert, sobol_map["vdw_BE"][0]] * DELTA_E0_MAX_J_MOL_VDW)/9.6e4
E_0_n = float(DELTA_E0_MAX_J_MOL - 2.0 * x_sobol[pert,sobol_map["N_BE"][0]] * DELTA_E0_MAX_J_MOL)/9.6e4

print(f"C: {E_0_c}\nO: {E_0_o}\nH: {E_0_h}\nVdw: {E_0_vdw}\nN: {E_0_n}")





