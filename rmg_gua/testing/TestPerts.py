# ### black box testing 
# check output from RMG run in a chemkin file and ensure perturbs were applied correctly
import rmgpy.chemkin as chmkn
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.data.surface import MetalDatabase
from torch.quasirandom import SobolEngine
import random
import pickle
import os
import yaml
import numpy as np


def test_thermo_pert(
    RMG_db_folder,
    output_path, 
    pert_yaml, 
    M,

):
    # perturbed model path 
    sobol_path = f"{output_path}sobol_map.pickle"
    with open(sobol_path, "rb") as f:
        perturb_values = pickle.load(f)
        
    # load the ranges
    sobol_range_path = f"{output_path}sobol_range_map.pickle"
    with open(sobol_range_path, "rb") as f:
        perturb_ranges = pickle.load(f)
    
    # config file with perturbation values
    with open(pert_yaml, 'r') as f:
        perts = yaml.safe_load(f)
  
    # Specify the path to the thermo library
    os.path.join(RMG_db_folder, "input", "thermo")
    thermo_library_path = os.path.join(RMG_db_folder, "input", "thermo")
    if not os.path.exists(thermo_library_path):
        raise OSError(f'Path to rules does not exist:\n{thermo_library_path}')

    # Specify the path to the metal binding energy library
    metal_path = os.path.join(RMG_db_folder,"input","surface")
    if not os.path.exists(metal_path):
        raise OSError(f'Path to rules does not exist:\n{metal_path}')

    # Thermo Perturbations
    # BEs that are in the metal database are in eV. vdw (in thermo db) is in
    # J/mol
    be_dict = {
        'C': perts["DELTA_E0_C"],
        'O': perts["DELTA_E0_O"],
        'H': perts["DELTA_E0_H"],
        'N': perts["DELTA_E0_N"],
        'Vdw': perts["DELTA_E0_VDW"],
    }

    ## sobol sequence testing
    # Create the pseudo randoms
    sobol = SobolEngine(dimension=300, scramble=True, seed=100)
    x_sobol = sobol.draw(M)
    
    # test that the ranges are correct
    for atom, be_range in be_dict.items():
        sobol_range = (perturb_ranges[atom][2] - perturb_ranges[atom][1])/2
        print(sobol_range,be_range)
        assert np.isclose(sobol_range,be_range)
        
        # next, get the range spanned by of all of the perturbed values
        # this could be difficult to test, make it fairly lenient for a/r tol
        sobol_val_range = (max(perturb_values[atom][1]) - min(perturb_values[atom][1]))/2
        print(sobol_val_range,be_range)
        assert np.isclose(sobol_val_range,be_range,rtol=1e-03, atol=1e-03)
    
    # test that the model is actually being perturbed
    
if __name__ == "__main__":
    
    # specify inputs 
    M = 5000 # number of runs of RMG to do
    N = 500 # number of slurm jobs to run at once
    unc_folder = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/"
    conda_path = os.path.join(unc_folder, "conda")
    RMG_base_folder = os.path.join(unc_folder, "RMG-Py")
    RMG_db_folder = os.path.join(unc_folder,"RMG-database")
    expt_yaml_file = os.path.join(unc_folder, "rmg_gua", "cantera" "all_experiments_reorg.yaml")

    rmg_unc_scripts_folder = os.path.join(unc_folder,"rmg_gua","rmg")
    output_path = "/scratch/blais.ch/methanol_results_2022_08_25_just_thermo_5000/"
    pert_yaml = os.path.join(unc_folder, "rmg_gua", "example", "perturb_groups.yaml")
    reference_input = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/rmg/input.py"
    
    
    test_thermo_pert(
        RMG_db_folder,
        output_path, 
        pert_yaml, 
        M,
    )