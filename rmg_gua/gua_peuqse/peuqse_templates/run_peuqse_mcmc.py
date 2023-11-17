import shutil
import os
import glob
import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput
import sys
repo_dir = __REPO_DIR__
sys.path.insert(0, repo_dir)
import rmg_gua.gua_peuqse.ct_simulation as ct_simulation
from rmg_gua.gua_peuqse.setup_peuqse import setup_userinput
project_path = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":

    print("running peuqse")
    # for continue sampling, peuqse doesn't reference the old pickle files right
    # until I can fix, moving manually for now. 
    try:
        files = glob.glob("./pickles/mpi_cached_files/*.pkl")

        for file in files:
            shutil.copy(file, "./pickles/")
        # continue_sampling if we have pickle files to copy
        if len(files) > 0:
            continue_sampling=True
        else: 
            continue_sampling=False
    except:
        continue_sampling=False
        print("could not copy files")
    
    print("running job")
    # setup our ct_simulation function
    ct_simulation.sim_init(project_path)

    # reduce space by excluding all of the rule alphas
    UserInput = setup_userinput(project_path, reduce_space="alpha")

    # saves all of the post burn in responses for ESS/mcmc. overwrites the old responses if you continue sampling. 
    UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] = True

    UserInput.parameter_estimation_settings['mcmc_length'] = 100
    # UserInput.parameter_estimation_settings['mcmc_nwalkers'] = 52
    UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] = __PARALLEL__
    UserInput.parameter_estimation_settings['mcmc_continueSampling'] = continue_sampling #IMPORTANT: for mcmc_parallel_sampling, you **must** put the continue in the UserInput and **before** the PE_object is created. Putting it as an argument in doMetropolisHastings won't work since all the previous logs will get deleted when the PE_object is initialized.

    UserInput.parameter_estimation_settings["mcmc_walkerInitialDistributionSpread"] = 0.25
    
    # trying filtering. auto is 2, so set it to 1. 
    UserInput.parameter_estimation_settings['mcmc_threshold_filter_coefficient'] = 1.0
    # UserInput.parameter_estimation_settings['mcmc_relative_step_length']
    
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)

    # mcmc_output = PE_object.doMetropolisHastings()
    mcmc_output = PE_object.doEnsembleSliceSampling()
    PEUQSE.save_PE_object(PE_object, "pe_object.dill")

    # PE_object = PEUQSE.load_PE_object(“FileBeingSavedAs.dill”)