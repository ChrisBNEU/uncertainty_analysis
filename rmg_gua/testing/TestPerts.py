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
import logging
import numpy as np
from tqdm import tqdm 



def get_diff(
    spec, 
    pert_spec, 
    mdb_dict, 
    pert_dict, 
    atom, 
    model_num, 
    eV = True,
    scale = 1,
    
):
    """
    simple function to check thermo
    spec = species from original model
    pert_spec = species from perturbed model
    mdb = metal database (for comparing original binding energy)
    pert_dict = sobol map of perturbed values
    atom = bound atom for linear scaling (e.g. O, H, or C)
    model_num = model number (e.g. 0-5000)
    eV = if true, convert enthalpies to eV
    scale = linear scaling multiplier for comparison (i.e. if O=C=X, 
    """

    # get the intended perturbed value. for vdw, original BE value is 0 
    # (no vdw BE in database, so arbitrarily reference at 0)
    if atom == 'Vdw':
        sobol_pert = -pert_dict[atom][1][model_num]
    else:
        sobol_pert = mdb_dict[atom].value - pert_dict[atom][1][model_num]
    
    # scale if comparison species is not atomic (i.e. c=x)
    sobol_pert = sobol_pert*scale
    
    # get the real pert value. if in eV, convert to eV
    if eV:
        convert = 9.6e4
    else: 
        convert = 1
        
    H_orig = spec.thermo.get_enthalpy(298)/convert
    H_pert = pert_spec.thermo.get_enthalpy(298)/convert
    actual_pert = H_orig - H_pert

    return sobol_pert, actual_pert

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
    
    mdb = MetalDatabase()
    mdb.load(metal_path)
    
    metal_vals_orig = mdb.libraries["surface"].entries[perts["metal"]].binding_energies
    # metal_vals_orig['C'].value
    # metal_vals_orig['O'].value, DELTA_E0_O, "BE_eV")
    # sobol_range_map["H"] = get_range(metal_vals_orig['H'].value, DELTA_E0_H, "BE_eV")
    # sobol_range_map["N"] = get_range(metal_vals_orig['N'].value, DELTA_E0_N, "BE_eV")
    # sobol_range_map["Vdw"] = get_range(0, DELTA_E0_VDW, "BE_J")

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
    
    # load original model
    orig_path = "/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/rmg_gua/baseline/"
    orig_chmkn_path = os.path.join(orig_path, "chemkin")
    orig_file_surf = os.path.join(orig_chmkn_path,"chem_annotated-surface.inp")
    orig_file_dict = os.path.join(orig_chmkn_path,"species_dictionary.txt")
    orig_species, orig_rxn = chmkn.load_chemkin_file(orig_file_surf, dictionary_path=orig_file_dict)
    
    # test that the model is actually being perturbed for all 0-5000
    for i in tqdm(range(0,M)):
        # load model
        model_str = str(i).zfill(4)
        chmkn_path = f"{output_path}run_{model_str}/chemkin/"
        file_surf = chmkn_path + "chem_annotated-surface.inp"
        file_gas = chmkn_path +  "chem_annotated-gas.inp"
        file_dict = chmkn_path +  "species_dictionary.txt"
        if os.path.exists(chmkn_path):
            
            pert_species, pert_rxn = chmkn.load_chemkin_file(file_surf, dictionary_path=file_dict)
            
            for num, spec_num in enumerate(orig_species):
                if spec_num.smiles == 'O=[Pt]':
                    oxygen_spec = spec_num
                    oxygen_pert_spec = pert_species[num]
                    sobol_pert_o, actual_pert_o = get_diff(
                        oxygen_spec, 
                        oxygen_pert_spec, 
                        metal_vals_orig, 
                        perturb_values, 
                        'O',
                        i,
                    )

                    if not np.isclose(sobol_pert_o,actual_pert_o, rtol=1e-04, atol=1e-02):
                        logging.warning(f"model {model_str} has bad_value for O perturbation ")
                        logging.warning(f"Perts (expected, actual): {sobol_pert_o} {actual_pert_o}")
                    # assert np.isclose(sobol_pert_o,actual_pert_o, rtol=1e-04, atol=1e-02)
                    # print("Perts (expected, actual): ", sobol_pert_o, actual_pert_o)
                    
                if spec_num.smiles == '[Pt]' and len(spec_num.molecule[0].atoms) > 1:
                    hydrogen_spec = spec_num
                    pert_hydrogen_spec = pert_species[num]
                    sobol_pert_h, actual_pert_h = get_diff(
                        hydrogen_spec, 
                        pert_hydrogen_spec, 
                        metal_vals_orig, 
                        perturb_values, 
                        'H',
                        i,
                    )
                    if not np.isclose(sobol_pert_h,actual_pert_h, rtol=1e-04, atol=1e-02):
                        logging.warning(f"model {model_str} has bad_value for H perturbation ")
                        logging.warning(f"Perts (expected, actual): {sobol_pert_h} {actual_pert_h}")
                    
                if spec_num.smiles == 'O=C=[Pt]':
                    # scaled by 0.5 because we are dealing with C=x, not atomic C
                    carbon_spec = spec_num
                    pert_carbon_spec = pert_species[num]
                    sobol_pert_c, actual_pert_c = get_diff(
                        carbon_spec , 
                        pert_carbon_spec, 
                        metal_vals_orig, 
                        perturb_values, 
                        'C',
                        i,
                        scale = 0.5
                    )
                    if not np.isclose(sobol_pert_c,actual_pert_c, rtol=1e-04, atol=1e-02):
                        logging.warning(f"model {model_str} has bad_value for C perturbation ")
                        logging.warning(f"Perts (expected, actual): {sobol_pert_c} {actual_pert_c}")
                    
                if spec_num.smiles == 'O=C=O.[Pt]':
                    vdw_spec = spec_num
                    pert_vdw_spec = pert_species[num]
                    sobol_pert_vdw, actual_pert_vdw = get_diff(
                        vdw_spec , 
                        pert_vdw_spec, 
                        metal_vals_orig, 
                        perturb_values, 
                        'Vdw',
                        i,
                        eV = False,
                    )
                    if not np.isclose(sobol_pert_vdw,actual_pert_vdw, rtol=1e-02, atol=1e-02):
                        logging.warning(f"model {model_str} has bad_value for vdw perturbation ")
                        logging.warning(f"Perts (expected, actual): {sobol_pert_vdw} {actual_pert_vdw}")
        else: 
            logging.warning(f"skipping run {model_str}")
                
        
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