import copy
import pickle
import shutil
import os
import math
import yaml
import pandas as pd
import time
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms

import numpy as np
import scipy.stats as stats
import sys
prefix = '/work/westgroup/ChrisB/_01_MeOH_repos/uncertainty_analysis/'
sys.path.append(prefix)
from gua_cantera.Spinning_basket_reactor.sbr import MinSBR
from rmg_gua.gua_peuqse.runtime_utilities import get_all_param_lists
from rmg_gua.gua_peuqse import ct_simulation
from rmg_gua.gua_peuqse.runtime_utilities import get_all_param_lists, setup_userinput

import PEUQSE
import cantera as ct

# conversion eV to j/Kmol
ev2jkm = 96.4869*1e6

class cpox_plots(): 

    def __init__(self, project_path,alternate_paths=None):

        self.load_data(project_path, alternate_paths=alternate_paths)
        
    def load_data(self,project_path, folder_name, by_species=False):
        """
        load all data from lsr ocp and dft runs
        """
        project_path_dft = os.path.join(project_path, folder_name)

        self.results_path = os.path.join(project_path, "config")
        self.params = get_all_param_lists(results_path=self.results_path, by_species=by_species)
        self.sim_wrappers = sim_wrapper
        self.sim_wrapper.results_path = self.results_path

        # these should all be the same, so load it only from one
        with open(os.path.join(self.results_path, "expt_vals.yaml"), "r") as f: 
            self.expt_vals = yaml.safe_load(f)

        # load uncertainties
        with open(os.path.join(self.results_path, "expt_unc.yaml"), "r") as f: 
            self.uncerts = yaml.safe_load(f)

        self.x_vals = self.expt_vals.pop("distance_mm")

        self.label_list = []
        self.value_list = []
        for key, value in self.expt_vals.items():
            self.label_list.append(key)
            self.value_list.append(value)
        

        # if not os.path.exists(os.path.join(project_path, "pickles/mcmc_post_burn_in_statistics.pkl")):
        pe_object = PEUQSE.load_PE_object(os.path.join(project_path, "pe_object.dill"))

        self.parameters_map = pe_object.map_parameter_set
        self.parameters_muap = pe_object.mu_AP_parameter_set
        self.parameters_stdap = pe_object.stdap_parameter_set
        
        # check dim of std prior. if 2-d, take diagonal to get the variance elements
        if len(pe_object.UserInput.std_prior.shape) ==2: 
            self.parameters_stdpriors[key] = np.diag(pe_object.UserInput.std_prior)
        else: 
            self.parameters_stdpriors[key] = pe_object.UserInput.std_prior
            
        with open(os.path.join(project_path, "pickles/mcmc_post_burn_in_statistics.pkl"), "rb") as input_file:
            burn_in_stats = pickle.load(input_file)
            print(key, ": ", len(burn_in_stats[5]))


        # simulate both the most probable and best fit set
        self.sim_wrappers[key].load_data(self.results_paths,by_spec=True,unc_src=unc_src, tprof=True, facet=211)
        self.output_maps[key] = self.sim_wrappers[key].monolith_wrapper(self.parameters_maps[key])

        self.sim_wrappers[key].load_data(self.results_paths[key],by_spec=True,unc_src=unc_src, tprof=True, facet=211)
        self.output_mu_aps[key] = self.sim_wrappers[key].monolith_wrapper(self.parameters_muaps[key])
        
        self.sim_wrappers[key].load_data(self.results_paths[key],by_spec=True,unc_src=unc_src, tprof=True, facet=211)
        self.output_initials[key] = self.sim_wrappers[key].monolith_wrapper([0]*13)



        # load the initial binding energies
        dft_yaml = os.path.join(repo_dir, "rh_dft_211", "rh_cpox_dft_211.yaml")
        lsr_yaml = os.path.join(repo_dir, "rh_211", "rh_cpox_lsr_211.yaml")
        ocp_yaml = os.path.join(repo_dir, "rh_ocp_211", "rh_cpox_ocp_211.yaml")

        self.yaml_paths_list = [dft_yaml, lsr_yaml, ocp_yaml]
        self.yaml_paths = dict(zip(self.project_names, self.yaml_paths_list))

        self.ct_sltns = {}
        for src, path in self.yaml_paths.items():
            gas = ct.Solution(path, 'gas')
            surf = ct.Interface(path, 'surface1', [gas])
            
            self.ct_sltns[src] = {"gas": gas, "surf": surf}


        # load all of the info needed for the observables (ch4 conversion, syngas selectivity, etc.)
        with open(os.path.join(self.project_paths[src], "config", "expt_vals.yaml"), "r") as f: 
            self.expt_vals = yaml.safe_load(f)


        self.distances = self.expt_vals.pop("distance_mm")
        self.temp_prof = self.expt_vals.pop("temp_profile")
        self.output_labels = list(self.expt_vals.keys())


        # we translate the labels to the species names so the script can interpret them
        self.label_lookup = {
            'ch4_profile':'CH4(2)',
            'co2_profile':'CO2(4)',
            'co_profile':'CO(7)',
            'h2_profile':'H2(6)',
            'h2o_profile':'H2O(5)',
            'o2_profile':'O2(3)',
        }
        self.species_labels = []
        for value in self.output_labels: 
            self.species_labels.append(self.label_lookup[value])

        # get the species names from the mechanism and make plot labels
        self.species_list = self.param_dicts["DFT"]["label_list"]
        self.species_list_pretty = [
            'H*',
            'CO$_2$*',
            'CO*',
            'CH$_4$*',
            'O* (211)',
            'CH$_2$*',
            'CH$_3$*',
            'CH*',
            'C*',
            'H$_2$*',
            'OH* (211)',
            'H$_2$O*',
            'CHO*'
        ]
        self.species_dict = dict(zip(self.species_list, self.species_list_pretty))

        # make nicer labels for plotting
        self.label_pretty = {}
        for label in self.label_list:
            if label == "ch4_profile": 
                self.label_pretty[label] = "$CH_4$ Profile" 
            if label == "co2_profile": 
                self.label_pretty[label] = "$CO_2$ Profile" 
            if label == "co_profile": 
                self.label_pretty[label] = "$CO$ Profile" 
            if label == "h2_profile": 
                self.label_pretty[label] = "$H_2$ Profile" 
            if label == "h2o_profile": 
                self.label_pretty[label] = "$H_{2}O$ Profile" 
            if label == "o2_profile": 
                self.label_pretty[label] = "$O_2$ Profile" 
            
        # identify any indices from the "simulation.calculate" function that are 
        # associated with temperature and distance. we will exclude those for now.
        self.observable_labels_sens = [
            "reference_syngas_selectivity", "reference_syngas_yield", "reference_co_sel", "reference_co_yield", 
            "reference_h2_sel", "reference_h2_yield", "reference_ch4_conv", "reference_full_oxidation_selectivity", 
            "reference_full_oxidation_yield", "reference_exit_temp",
            "reference_peak_temp", "reference_peak_temp_dist", "reference_o2_conv", "reference_max_ch4_conv", "reference_dist_to_50_ch4_conv"]

        self.remove_indices = []
        self.keep_obs_labels = []
        for index, label in enumerate(self.observable_labels_sens): 
            if "temp" in label or "dist" in label:
                self.remove_indices.append(index)
            else: 
                self.keep_obs_labels.append(label.replace("reference_", ""))
    


    def make_yamls(self, map_yaml_path):
        """ 
        generate output yaml files

        yaml_path: path to save the yaml files
        """
        # generate the yaml files for the map dataset
        mech_yaml = os.path.join(repo_dir, "rh_dft_211", "rh_cpox_dft_211.yaml")


        self.yaml_paths = dict(zip(self.project_names, yaml_path_list))

        yaml_path = self.yaml_paths[src]
        project_path = self.project_paths[src]
        results_path = os.path.join(project_path, "config")
        # phase definition with 'density' for the surface phase.
        # this throws an error in 3.0, so remove
        with open(yaml_path, "r") as f: 
            mech_yaml = yaml.load(f, Loader=yaml.FullLoader)

        for phnum, phase in enumerate(mech_yaml["phases"]): 
            if phase["name"] == "surface1" and "density" in phase["state"].keys():
                mech_yaml["phases"][phnum]["state"].pop("density")

        new_yaml_path = "./new_yaml.yaml"
        with open(new_yaml_path, "w") as f: 
            yaml.dump(mech_yaml,f)

        # load perturbations for mechanism
        # post burn in statistics has: 
        # self.map_parameter_set, 
        # self.mu_AP_parameter_set, 
        # self.stdap_parameter_set, 
        # self.evidence, 
        # self.info_gain, 
        # self.post_burn_in_samples, 
        # self.post_burn_in_log_posteriors_un_normed_vec
        with open(os.path.join(project_path, "pickles/mcmc_post_burn_in_statistics.pkl"), "rb") as input_file:
            burn_in_stats = pickle.load(input_file)

        parameters_map = burn_in_stats[0]
        # parameters_muap = burn_in_stats[1]
        # parameters_stdap = burn_in_stats[2]

        perts = parameters_map
        modified_mech = MinSBR.setup_ct_solution(new_yaml_path, pert=perts, results_path=results_path, by_species=True, debug=False)

        yml_write = ct.YamlWriter()
        yml_write.set_header(modified_mech["gas"])
        yml_write.add_solution(modified_mech["gas"])
        yml_write.add_solution(modified_mech["surf"])
        yml_write.to_file(os.path.join(map_yaml_path, f"./{src}_map_mechanism.yaml"))

        # remove the intermediate yaml file
        os.remove(new_yaml_path)

#     def get_binding_energies(self, use_zpe_be_start=False, facet="211"):
#         """
#         get binding energies for species

#         use_zpe_be_start: if true, use the zpe corrected binding energies for the starting point
#         else, use the h298 from the cantera file
#         """
#         # load zpe corrected binding energy
#         if use_zpe_be_start:
#             zpe_correction_dft = os.path.join(repo_dir, "plotting", "scaled_bes", f"dft_be_dict_{facet}.yaml")
#             zpe_correction_lsr = os.path.join(repo_dir, "plotting", "scaled_bes", f"lsr_be_dict_{facet}.yaml")
#             zpe_correction_ocp = os.path.join(repo_dir, "plotting", "scaled_bes", f"ocp_be_dict_{facet}.yaml")
#             with open(zpe_correction_dft, "r") as f: 
#                 zpe_correction_dft = yaml.safe_load(f)
#             with open(zpe_correction_lsr, "r") as f:
#                 zpe_correction_lsr = yaml.safe_load(f)
#             with open(zpe_correction_ocp, "r") as f:
#                 zpe_correction_ocp = yaml.safe_load(f)

#             zpe_corrections = {"DFT": zpe_correction_dft, "LSR": zpe_correction_lsr, "ML": zpe_correction_ocp}

#         # get the gas phase precursor hf298s loaded
#         with open(os.path.join(repo_dir, "rh", "gas_precursor_hf298.yaml"), "r") as f: 
#             gas_prec_hf298s = yaml.safe_load(f)

#         self.spec_start_h298s = {}
#         self.spec_uncs = {}
#         self.spec_uncs_ub = {}
#         self.spec_uncs_lb = {}
#         self.params_map_scaled = {}
#         self.params_muap_scaled = {}
#         self.params_stdap_scaled = {}

#         for src, sltn in self.ct_sltns.items():
#             spec_start_h298 = {}
#             map_scaled = {}
#             muap_scaled = {}
#             stdap_scaled = {}
            
#             map_dict = dict(zip(self.param_dicts[src]["label_list"], self.parameters_maps[src]))
#             muap_dict = dict(zip(self.param_dicts[src]["label_list"], self.parameters_muaps[src]))
#             stdap_dict = dict(zip(self.param_dicts[src]["label_list"], self.parameters_stdaps[src]))
            
#             with open(os.path.join(self.results_paths[src], "be_unc.yaml")) as f: 
#                 self.spec_uncs[src] = yaml.safe_load(f)
                
#             self.spec_uncs_ub[src] = copy.deepcopy(self.spec_uncs[src])
#             self.spec_uncs_lb[src] = copy.deepcopy(self.spec_uncs[src])
                    
#             for param in self.spec_uncs[src].keys():
#                 # scale out uncertainties to eV 
#                 self.spec_uncs[src][param] = self.spec_uncs[src][param]/ev2jkm
                
#                 # get value for initial h298 value
#                 spec = self.ct_sltns[src]["surf"].species(param)

#                 if use_zpe_be_start: 
#                     spec_start_h298[param] = zpe_corrections[src][param]
#                 else:
#                     spec_start_h298[param] = spec.thermo.h(298)/(ev2jkm) - gas_prec_hf298s[param]
                
#                 # get upper and lower bound for uncertainty
#                 self.spec_uncs_ub[src][param] = spec_start_h298[param]+self.spec_uncs[src][param]
#                 self.spec_uncs_lb[src][param] = spec_start_h298[param]-self.spec_uncs[src][param]
                
#                 # scale the map, muap, and std ap using the starting points specified 
#                 map_scaled[param] = map_dict[param]/(ev2jkm) + spec_start_h298[param] 
#                 muap_scaled[param] = muap_dict[param]/(ev2jkm) + spec_start_h298[param] 
#                 stdap_scaled[param] = stdap_dict[param]/(ev2jkm) + spec_start_h298[param] 
            
#             self.spec_start_h298s[src] = spec_start_h298
#             self.params_map_scaled[src] = map_scaled
#             self.params_muap_scaled[src] = muap_scaled
#             self.params_stdap_scaled[src] = stdap_scaled

    
#     def load_gsa_data(self,gsa_data_path=None, no_exp_err=False):


#         def conflicts(a, b):
#             """ 
#             helper function for picking out the post burn in data from the aggregate data
#             """
#             from collections import defaultdict
#             elt2ix = defaultdict(list)
#             for i, elt in enumerate(a):
#                 elt2ix[elt].append(i)
#             for j, elt in enumerate(b):
#                 if elt in elt2ix:
#                     for i in elt2ix[elt]:
#                         yield i, j
#         # load the gsa output file: 
#         self.gsa_outputs = {}
#         self.gsa_outputs_logps = {}
#         self.gsa_outputs_bes = {}
#         self.gsa_outputs_res = {}
#         for src, path in self.project_paths.items(): 
#             with open(os.path.join(path, "pickles", "gsa_output_all.pkl"), "rb") as f:
#                 self.gsa_outputs[src] = pickle.load(f)
            
#             gsa_output_np = np.array(self.gsa_outputs[src], dtype=object)
            
#             gsa_output_np_copy = copy.deepcopy(gsa_output_np)
            
#             # just give everything that is less than 0 or NAN to 0 probability
#             for idx, data in enumerate(gsa_output_np_copy): 
#                 if not isinstance(data[0], np.float64) and len(data[0]) == 6:
#                     if  data[0][0] == np.NAN or data[0][0] == np.NINF or data[0][0] < 0: 
#                         gsa_output_np[idx, 0] = 0.0
#                 elif isinstance(data[0], np.float64):
#                     if data[0] == np.NAN or data[0] == np.NINF or data[0] < 0: 
#                         gsa_output_np[idx, 0] = 0.0
#                 else: 
#                     raise ValueError("Data is not the right shape: ", data[0].shape)
                        
#             print(src)  
#             self.gsa_outputs_logps[src] = np.array(list(gsa_output_np[:,0]), dtype=np.ndarray)
#             self.gsa_outputs_bes[src] = np.array(list(gsa_output_np[:,1]), dtype=np.ndarray)
#             self.gsa_outputs_res[src] = np.array(list(gsa_output_np[:,2]), dtype=np.ndarray)
        
#         self.gsa_outputs_logps_sanit = {}
#         self.burn_in_outputs_np = {}
#         self.burn_in_outputs_std = {}
#         self.burn_in_bes = {}
#         self.burn_in_obs = {}
#         self.burn_in_obs_notemp = {}

#         # loop through projects and 
#         for src in self.project_names: 
#             project_path = self.project_paths[src]
#             pe_object = PEUQSE.load_PE_object(os.path.join(project_path, "pe_object.dill"))
            
#             # get all unique binding energy inputs. then, get the matching output values
#             gsa_bes_contig = np.ascontiguousarray(self.gsa_outputs_bes[src], dtype=np.float64)
#             gsa_bes_sanit, uniq_bes_indx = np.unique(gsa_bes_contig,axis=0,return_index=True)
#             self.gsa_outputs_sanit = self.gsa_outputs_res[src][uniq_bes_indx]
            
#             # drop any output rows that are "None", and the corresponding inputs
#             keeprows = [i for i,v in enumerate(self.gsa_outputs_sanit) if not isinstance(v, type(None))]
#             gsa_bes_sanit = gsa_bes_sanit[keeprows]
#             self.gsa_outputs_sanit = self.gsa_outputs_sanit[keeprows]

#             # do the same for logps
#             self.gsa_outputs_logps_sanit[src] = self.gsa_outputs_logps[src][keeprows]
            
#             assert len(self.gsa_outputs_sanit) == len(gsa_bes_sanit) 
#             assert len(self.gsa_outputs_sanit) == len(self.gsa_outputs_logps_sanit[src])
            
#             self.burn_in_bes[src] = pe_object.post_burn_in_samples
            
#             list1 = [tuple(x) for x in gsa_bes_sanit.astype(float).tolist()]
#             list2 = [tuple(x) for x in pe_object.post_burn_in_samples.astype(float).tolist()]
            
#             # overlap = conflicts(list1, list2)
#             burn_in_outputs = []
#             self.burn_in_obs_list = []
#             conflict_indices = []
#             conflict_gen = conflicts(list1, list2)

#             # get original shape of data for one parameter set
#             # need a work around for no expt error data. 
#             if no_exp_err:
#                 x,y = 6,22
#             else: 
#                 x,y = self.gsa_outputs_sanit[0].shape

#             burn_in_outputs_list = []
#             for h, (i,j) in enumerate(conflict_gen): 
#                 # if len(burn_in_outputs) == 0: 
#                 # # conflict_indices.append([i,j])
#                 burnin_output = self.gsa_outputs_sanit[i].reshape(1,x,y)
#                 burn_in_outputs.append(burnin_output)
                
#                 # get calculated outputs such as ch4 conversion, selectivity, etc.
#                 observables = simulation.calculate([burnin_output[0], self.species_labels, self.distances, self.distances, None])
#                 self.burn_in_obs_list.append(np.array(observables).reshape(1, len(observables)))
            
#             # string list of burn in simulated outputs to a 3-d np array
#             self.burn_in_outputs_np[src] = np.concatenate(burn_in_outputs, axis=0)
#             self.burn_in_obs[src] = np.concatenate(self.burn_in_obs_list)
#             self.burn_in_obs_notemp[src] = np.delete(self.burn_in_obs[src], self.remove_indices, axis=1)
            
            
#             # indexing: x = burn in run, y = trend (ch4, o2, etc), z = point in reactor (spatially)
#             # construct a y*z dim array for the standard deviations 
#             t1 = time.time()
#             self.burn_in_outputs_std[src] = np.zeros((self.burn_in_outputs_np[src].shape[1], self.burn_in_outputs_np[src].shape[2]))
#             for z in range(self.burn_in_outputs_np[src].shape[2]): 
#                 for y in range(self.burn_in_outputs_np[src].shape[1]):
#                     if np.isnan(np.sum(self.burn_in_outputs_np[src][:,y,z])):
#                         self.burn_in_outputs_std[src][y,z] = np.nanstd(self.burn_in_outputs_np[src][:,y,z].astype(float))
#                     else:
#                         self.burn_in_outputs_std[src][y,z] = self.burn_in_outputs_np[src][:,y,z].std()
#             t2 = time.time()
#             print("std calc time: ", t2-t1)

#     def plot_outputs(self, plot_path=None, plot_name = None, no_exp_dat_path=None):
#         """ 
#         plot the prior and posterior distributions for the binding energies, 
#         along with the map values and the starting point values for be

#         plot path: folder to save the plots in

#         returns the figure and axis objects for use in jupyter
#         """

#         # plot most probable and best fit set
#         plt.rcParams["figure.figsize"] = [15,3]
#         # fig, ax = plt.subplots(len(self.value_list)-1,3)
        
#         y_axis_labels = {
#             'ch4_profile': 'CH$_4$',
#             'co2_profile': 'CO$_2$',
#             'co_profile': 'CO',
#             'h2_profile': 'H$_2$',
#             'h2o_profile': 'H$_2$O',
#             'o2_profile': 'O$_2$'
#             }


#         plts = []
#         for axis, op in enumerate(self.value_list):
#             if self.label_list[axis] != "temp_profile":
#                 fig, ax = plt.subplots(1,3)
                
#                 # loop over lsr, dft, ocp results
#                 for idx, src in enumerate(self.project_names):
#                     # ax[idx].set_title(f"{src} " + self.label_pretty[self.label_list[axis]])
#                     ax[idx].scatter(self.x_vals, self.output_initials[src][axis], label=f"Initial", color="g")
#                     ax[idx].scatter(self.x_vals, self.output_maps[src][axis], label=f"MAP", color="orange")
#                     # ax[axis,idx].scatter(self.x_vals, self.output_mu_aps[src][axis], label=f"{src} - mu_ap", color="orange")
#                     # ax[axis,idx].scatter(self.x_vals, self.output_std_aps[src][axis], label=f"{src} - std_ap")
#                     ax[idx].plot(self.x_vals, op, label="Experiment", linestyle="-", linewidth=1.5, marker="o")
#                     # ax[axis,idx].errorbar(self.x_vals, op, yerr=self.uncerts[self.label_list[axis]], capsize=2.5, capthick=1)
#                     # ax[axis].errorbar(self.x_vals, op, yerr=flat_self.uncerts[self.label_list[axis]], capsize=2.5, capthick=1)
#                     # ax[axis].errorbar(self.x_vals, op, yerr=self.uncerts[self.label_list[axis]], capsize=2.5, capthick=1)
#                     # ax[axis,idx].errorbar(self.x_vals, self.output_mu_aps[src][axis], yerr=self.burn_in_outputs_std[axis], capsize=2.5, capthick=1)

#                     # trim lower error bars if they go below 0
#                     # *2 for 2 sigma (95% ci)
#                     exp_upper_error = np.array(op) + np.array(self.uncerts[self.label_list[axis]])*2
#                     exp_lower_error = np.array(op) - np.array(self.uncerts[self.label_list[axis]])*2
#                     exp_lower_error[exp_lower_error < 0] = 0

#                     sim_upper_error = self.output_maps[src][axis] + self.burn_in_outputs_std[src][axis]*2
#                     sim_lower_error = self.output_maps[src][axis] - self.burn_in_outputs_std[src][axis]*2
#                     sim_lower_error[sim_lower_error < 0] = 0

#                     ax[idx].fill_between(
#                         self.x_vals, 
#                         # self.output_mu_aps[src][axis] - self.burn_in_outputs_std[axis],
#                         # self.output_mu_aps[src][axis] + self.burn_in_outputs_std[axis],
#                         sim_lower_error,
#                         sim_upper_error,
#                         color="sandybrown",
#                         alpha= 0.3,
#                         label="2$\sigma$ MAP",
#                         )
#                     ax[idx].fill_between(
#                         self.x_vals, 
#                         exp_lower_error,
#                         exp_upper_error,
#                         color="powderblue",
#                         alpha= 0.3,
#                         label="2$\sigma$ Experiment"
#                         )
                    
#                     # if we specify it, include uncertainty for initial data when we exclude experimental data uncertainty
#                     if no_exp_dat_path: 
#                         with open(no_exp_dat_path, "rb") as f:
#                             initial_point_err = pickle.load(f)
                            
#                         init_upper_error = self.output_initials[src][axis] + initial_point_err[src][axis]*2
#                         init_lower_error = self.output_initials[src][axis] - initial_point_err[src][axis]*2
#                         init_lower_error[init_lower_error < 0] = 0
#                         ax[idx].fill_between(
#                             self.x_vals, 
#                             init_upper_error,
#                             init_lower_error,
#                             color="g",
#                             alpha= 0.1,
#                             label="2$\sigma$ Initial"
#                             )
                        
                    
#                     if idx == 0:
#                         ax[idx].set_ylabel(f"{y_axis_labels[self.label_list[axis]]} flow mol/min")
                    
#                     ax[idx].set_xlabel("Reactor coordinate (mm)")
                    
#                     if idx ==2: 
#                         ax[idx].legend(bbox_to_anchor=(1.5,0.9))
#                         buff = 0.15
#                         # at end, readjust the y limits to make sure they are all identical
#                         ymax = max([ax[0].get_ylim()[1], ax[1].get_ylim()[1], ax[2].get_ylim()[1]])
#                         ymin = min([ax[0].get_ylim()[0], ax[1].get_ylim()[0], ax[2].get_ylim()[0]])
#                         ax[0].set_ylim(ymin, ymax+ymax*buff)
#                         ax[1].set_ylim(ymin, ymax+ymax*buff)
#                         ax[2].set_ylim(ymin, ymax+ymax*buff)
#                     fig.canvas.draw()
                        
                
#                 if plot_path is not None: 
#                     save_path = os.path.join(plot_path, f"{self.label_list[axis]}.pdf")
#                     fig.savefig(save_path,bbox_inches="tight")
                    
#                 plts.append((fig, ax)) 

#         return plts   
    
#     ###########################################

#     def plot_o_oh_111_211(self, paper_plots_111, plot_path=None):
#         """
#         plot O and OH on 111 and 211
#         load the 
#         """

#         # load the 111 data if it is not already loaded
#         self.plots_111 = paper_plots_111

#         # figure size
#         plt.rcParams["figure.figsize"] = [15,3]

#         # set font size
#         plt.rcParams.update({'font.size': 12})

#         # list of fig ax objects to return after, for use in jupyter
#         plts = []

#         # grab posterior stdev values
#         stdaps_dict = {}
#         stdpriors_dict = {}
#         for src in self.project_names:
#             stdaps_dict[src] = dict(zip(self.param_dicts[src]["label_list"], self.parameters_stdaps[src]))
#             stdpriors_dict[src] = dict(zip(self.param_dicts[src]["label_list"], self.parameters_stdpriors[src]))

#         # outer = gridspec.GridSpec(nrows=2, ncols=2, wspace=0.2, hspace=0.5)


#         self.species_list_pretty = [
#                     'H*',
#                     'CO$_2$*',
#                     'CO*',
#                     'CH$_4$*',
#                     'O* (211)',
#                     'CH$_2$*',
#                     'CH$_3$*',
#                     'CH*',
#                     'C*',
#                     'H$_2$*',
#                     'OH* (211)',
#                     'H$_2$O*',
#                     'CHO*'
#                 ]

#         row = 0
#         for datidx, (species, pretty_species) in enumerate(self.species_dict.items()):

#             # just do ocp species
#             if species != "OX(25)" and species != "OHX(31)":
#                 continue
#             fig = plt.figure(figsize=(5, 3))
#             inner = gridspec.GridSpec(3, 1, wspace=0.1, hspace=0.1)                    
#             # inner = gridspec.GridSpecFromSubplotSpec(3, 1,
#             #     subplot_spec=outer[row,1], wspace=0.1, hspace=0.1)

#             xmaxs = []
#             xmins = []

#             # define the max and min for the xaxis 
#             for src in self.project_names:
#                 xmaxs.append(self.spec_uncs_ub[src][species])
#                 xmins.append(self.spec_uncs_lb[src][species])

#                 # also append the max values from map, muap
#                 xmaxs.append(self.params_map_scaled[src][species])
#                 xmins.append(self.params_map_scaled[src][species])

#                 xmaxs.append(self.params_muap_scaled[src][species])
#                 xmins.append(self.params_muap_scaled[src][species])

#                 # also need to look at the new uncertainties for 
#                 # upper and lower bounds for the plot
#                 xmaxs.append(self.params_map_scaled[src][species] + stdaps_dict[src][species]/ev2jkm)
#                 xmins.append(self.params_map_scaled[src][species] - stdaps_dict[src][species]/ev2jkm)

#             # make the buffer 2*sigma so we can see the distribution better
#             # so far ocp has the largest uncertainty, So just hardcoding that
#             buffer = stdaps_dict["ML"][species]/ev2jkm*2 


#         #     xmax = max(xmaxs) + buffer
#         #     xmin = min(xmins) - buffer
#             # hardcoding x limits for this plot: 
#             if species == "OX(25)": 
#                 xmin = -6.5
#                 xmax = -3.0
#             else: 
#                 xmin = -4.5
#                 xmax = -1.0

#             # plot dft, lsr, ocp as the columns
#             axes = []
#             for col, src in enumerate(self.project_names):
#                 initval = self.spec_start_h298s[src][species]
#                 mapval = self.params_map_scaled[src][species]
#                 muapval = self.params_muap_scaled[src][species] 

#                 # our prior uncertainties
#                 ub = self.spec_uncs_ub[src][species]
#                 lb = self.spec_uncs_lb[src][species]

#                 # the posterior uncertainties
#                 ubstdap = mapval + stdaps_dict[src][species]/ev2jkm
#                 lbstdap = mapval - stdaps_dict[src][species]/ev2jkm

#                 ax = plt.Subplot(fig, inner[col])
#                 axes.append(ax)
#                 show_lit_vals = True

#                 # ax.axvline(x=mapval, color='r', label = f"map")

#                 # # source: abild pedersen paper
#                 if species == "OHX(31)":
#                     plot_name = "211ohx"
#                     ylabel = "Rh(211)"
#                     # ax.axvline(x=-3.26, color='b', label = f"literature (211)")

#                 # source: abild pedersen paper
#                 elif species == "OX(25)":
#                     plot_name = "211ox"
#                     ylabel = ""
#                     # ax.axvline(x=-4.9, color='b', label = f"literature (211)")

#                 # make the xmin/xmax the max and min of the parameter uncertainty_vals
#                 ax.set_xlim(left=xmin, right=xmax)

#                 # plot normal distribution 
#                 mu = initval
#                 sigma = stdpriors_dict[src][species]/ev2jkm
#                 # ax.axvline(x=initval+sigma, color='b')
#                 # ax.axvline(x=initval-sigma, color='b')
#                 print("211 ", src, species, sigma)
#                 x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
#                 ax.fill_between(x, stats.norm.pdf(x, mu, sigma), alpha=0.2, facecolor="black", color=None, label="prior distribution" )

#                 # plot distribution of the post burn in data: 
#                 x_data = self.burn_in_bes[src][:,datidx]/ev2jkm + self.spec_start_h298s[src][species]
#                 ax.hist(x=x_data, bins = 100, density=True, alpha = 0.2, color= "m", label="posterior distribution",) 

#                 # create a secondary y axis 
#                 ax2 = ax.twinx()
#                 ax2.set_ylabel(src,rotation='horizontal')
#                 ax2.set_yticklabels([])
#                 ax2.set_yticks([])
#                 ax2.yaxis.set_label_coords(1.05, 0.7)

#                 # only show y-axis label if it is second row
#                 if col == 0:
#                     ax.get_xaxis().set_visible(False)
#                     ax.get_yaxis().set_visible(False)
#                 elif col ==1: 
#                     ax.get_xaxis().set_visible(False)
#                     ax.set_ylabel(ylabel, fontsize=16)
#                     ax.get_yaxis().set_visible(True)
#                 elif col == 2:
#                     ax.get_xaxis().set_visible(True)
#                     ax.set_xlabel("Binding Energy (eV)")
#                     ax.get_yaxis().set_visible(False)
#                     # ax.legend(bbox_to_anchor=(1.1, 3.5), loc="upper left")

#                 # hide the y axis
#                 ax.set_yticklabels([])
#                 ax.set_yticks([])
#                 fig.add_subplot(ax)

#             # save only the legend: 
#             if row == 0 and col == 1: 
#                 legend = ax.get_legend()
#                 bbox = mpl.transforms.Bbox([[497.08333333333337, -12.718749999999986], [732.7083333333334, 110.94791666666667]])
#                 fig.canvas.draw()
#                 bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

#                 # add a small buffer to the legend bounding box, so we don't get the legend outline
#                 bbox.x0 = bbox.x0 + abs(bbox.x0)*0.01
#                 bbox.x1 = bbox.x1 - abs(bbox.x1)*0.01
#                 bbox.y0 = bbox.y0 + abs(bbox.y0)*0.05
#                 bbox.y1 = bbox.y1 - abs(bbox.y1)*0.05
#                 # fig.savefig('legend.pdf', dpi="figure",bbox_inches=bbox)


#             if plot_path is not None: 
#                 save_path = os.path.join(plot_path, f"{plot_name}.pdf")
#                 fig.savefig(save_path, bbox_inches="tight")
#             row+=1

#         row=0

#         stdaps_dict = {}
#         stdpriors_dict = {}
#         for src in self.project_names:
#             stdaps_dict[src] = dict(zip(self.plots_111.param_dicts[src]["label_list"], self.plots_111.parameters_stdaps[src]))
#             stdpriors_dict[src] = dict(zip(self.plots_111.param_dicts[src]["label_list"], self.plots_111.parameters_stdpriors[src]))

#         for datidx, (species, pretty_species) in enumerate(self.plots_111.species_dict.items()):

#             # just do ocp species
#             if species != "OX(25)" and species != "OHX(31)":
#                 continue

#             inner = gridspec.GridSpec(3, 1, wspace=0.1, hspace=0.1)                    
#             # inner = gridspec.GridSpecFromSubplotSpec(3, 1,
#             #     subplot_spec=outer[row,0], wspace=0.1, hspace=0.1)

#             fig = plt.figure(figsize=(5, 3))
#             xmaxs = []
#             xmins = []

#             # define the max and min for the xaxis 
#             for src in self.plots_111.project_names:
#                 xmaxs.append(self.plots_111.spec_uncs_ub[src][species])
#                 xmins.append(self.plots_111.spec_uncs_lb[src][species])

#                 # also append the max values from map, muap
#                 xmaxs.append(self.plots_111.params_map_scaled[src][species])
#                 xmins.append(self.plots_111.params_map_scaled[src][species])

#                 xmaxs.append(self.plots_111.params_muap_scaled[src][species])
#                 xmins.append(self.plots_111.params_muap_scaled[src][species])

#                 # also need to look at the new uncertainties for 
#                 # upper and lower bounds for the plot
#                 xmaxs.append(self.plots_111.params_map_scaled[src][species] + stdaps_dict[src][species]/ev2jkm)
#                 xmins.append(self.plots_111.params_map_scaled[src][species] - stdaps_dict[src][species]/ev2jkm)

#             # make the buffer 2*sigma so we can see the distribution better
#             # so far ocp has the largest uncertainty, So just hardcoding that
#             buffer = stdaps_dict["ML"][species]/ev2jkm*2 


#             # xmax = max(xmaxs) + buffer
#             # xmin = min(xmins) - buffer
#             # hardcoding x limits for this plot: 
#             if species == "OX(25)": 
#                 xmin = -6.5
#                 xmax = -3.0
#             elif species == "OHX(31)": 
#                 xmin = -4.5
#                 xmax = -1.0

#             # plot dft, lsr, ocp as the columns
#             axes = []
#             for col, src in enumerate(self.plots_111.project_names):
#                 initval = self.plots_111.spec_start_h298s[src][species]
#                 mapval = self.plots_111.params_map_scaled[src][species]
#                 muapval = self.plots_111.params_muap_scaled[src][species] 

#                 # our prior uncertainties
#                 ub = self.plots_111.spec_uncs_ub[src][species]
#                 lb = self.plots_111.spec_uncs_lb[src][species]

#                 # the posterior uncertainties
#                 ubstdap = mapval + stdaps_dict[src][species]/ev2jkm
#                 lbstdap = mapval - stdaps_dict[src][species]/ev2jkm

#                 ax = plt.Subplot(fig, inner[col])
#                 axes.append(ax)
#                 show_lit_vals = True

#                 # ax.axvline(x=mapval, color='r', label = f"map")

#                 # # source: abild pedersen paper
#                 if species == "OHX(31)":
#                     # ylabel = "OH* (111)"
#                     ylabel = "Rh(111)"
#                     fig.suptitle('OH*', fontsize=16)
#                     plot_name = "111ohx"
#                     # ax.axvline(x=-2.39, color='b', label = f"literature (111)")


#                 # source: abild pedersen paper
#                 elif species == "OX(25)":
#                     ylabel = ""
#                     fig.suptitle('O*', fontsize=16)
#                     plot_name = "111ox"
#                     # ax.axvline(x=-4.66, color='b', label = f"literature (111)")

#                 # make the xmin/xmax the max and min of the parameter uncertainty_vals
#                 ax.set_xlim(left=xmin, right=xmax)

#                 # plot normal distribution 
#                 mu = initval
#                 # sigma = stdaps_dict[src][species]/ev2jkm
#                 sigma = stdpriors_dict[src][species]/ev2jkm
#                 # ax.axvline(x=initval+sigma, color='b')
#                 # ax.axvline(x=initval-sigma, color='b')

#                 x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
#                 ax.fill_between(x, stats.norm.pdf(x, mu, sigma), alpha=0.2, facecolor="black", color=None, label="prior distribution" )

#                 # plot distribution of the post burn in data: 
#                 x_data = self.plots_111.burn_in_bes[src][:,datidx]/ev2jkm + self.plots_111.spec_start_h298s[src][species]
#                 ax.hist(x=x_data, bins = 100, density=True, alpha = 0.2, color= "m", label="posterior distribution",) 

#                 # create a secondary y axis 
#                 ax2 = ax.twinx()
#                 ax2.set_ylabel(src,rotation='horizontal')
#                 ax2.set_yticklabels([])
#                 ax2.set_yticks([])
#                 ax2.yaxis.set_label_coords(1.05, 0.7)

#                 # only show y-axis label if it is second row
#                 if col == 0:
#                     ax.get_xaxis().set_visible(False)
#                     ax.get_yaxis().set_visible(False)
#                 elif col ==1: 
#                     ax.get_xaxis().set_visible(False)
#                     ax.set_ylabel(ylabel, fontsize = 16)
#                     ax.get_yaxis().set_visible(True)
#                 elif col == 2:
#                     ax.get_xaxis().set_visible(True)
#                     ax.set_xlabel("Binding Energy (eV)")
#                     ax.get_yaxis().set_visible(False)
#                     # ax.legend(bbox_to_anchor=(1.1, 3.5), loc="upper left")

#                 # hide the y axis
#                 ax.set_yticklabels([])
#                 ax.set_yticks([])
#                 fig.add_subplot(ax)

#             # save only the legend: 
#             if row == 0 and col == 1: 
#                 legend = ax.get_legend()
#                 bbox = mpl.transforms.Bbox([[497.08333333333337, -12.718749999999986], [732.7083333333334, 110.94791666666667]])
#                 fig.canvas.draw()
#                 bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

#                 # add a small buffer to the legend bounding box, so we don't get the legend outline
#                 bbox.x0 = bbox.x0 + abs(bbox.x0)*0.01
#                 bbox.x1 = bbox.x1 - abs(bbox.x1)*0.01
#                 bbox.y0 = bbox.y0 + abs(bbox.y0)*0.05
#                 bbox.y1 = bbox.y1 - abs(bbox.y1)*0.05

#             if plot_path is not None: 
#                 save_path = os.path.join(plot_path, f"{plot_name}.pdf")
#                 fig.savefig(save_path, bbox_inches="tight")
#             row+=1

    
# ###########################################
#     def plot_bes(self, plot_path=None, plot_name=None): 
#         """ 
#         plot the posterior distribution for the molar flow rate data
#         """
        
#         if not plot_name: 
#             plot_name = "figure"
#         # load the 
#         # figure size
#         plt.rcParams["figure.figsize"] = [15,3]

#         # set font size
#         plt.rcParams.update({'font.size': 12})

#         # list of fig ax objects to return after, for use in jupyter
#         plts = []

#         # grab posterior stdev values
#         stdaps_dict = {}
#         stdpriors_dict = {}
#         for src in self.project_names:
#             stdaps_dict[src] = dict(zip(self.param_dicts[src]["label_list"], self.parameters_stdaps[src]))
#             stdpriors_dict[src] = dict(zip(self.param_dicts[src]["label_list"], self.parameters_stdpriors[src]))
            
            
#         for row, (species, pretty_species) in enumerate(self.species_dict.items()):
#             fig = plt.figure(figsize=(5, 3))
#             inner = gridspec.GridSpec(3, 1, wspace=0.1, hspace=0.1)                    
            
#             xmaxs = []
#             xmins = []

#             # define the max and min for the xaxis 
#             for src in self.project_names:
#                 xmaxs.append(self.spec_uncs_ub[src][species])
#                 xmins.append(self.spec_uncs_lb[src][species])
                
#                 # also append the max values from map, muap
#                 xmaxs.append(self.params_map_scaled[src][species])
#                 xmins.append(self.params_map_scaled[src][species])
            
#                 xmaxs.append(self.params_muap_scaled[src][species])
#                 xmins.append(self.params_muap_scaled[src][species])
                
#                 # also need to look at the new uncertainties for 
#                 # upper and lower bounds for the plot
#                 xmaxs.append(self.params_map_scaled[src][species] + 2*stdpriors_dict[src][species]/ev2jkm)
#                 xmins.append(self.params_map_scaled[src][species] - 2*stdpriors_dict[src][species]/ev2jkm)
            
#             # make the buffer 2*sigma so we can see the distribution better
#             # so far ocp has the largest uncertainty, So just hardcoding that
#             buffer = 2*stdaps_dict["ML"][species]/ev2jkm*2 
            
#             # hardcoding max/min for oh and o
#             if species != "OHX(31)" and species !="OX(25)": 
#                 xmax = max(xmaxs) + buffer
#                 xmin = min(xmins) - buffer
                
#             elif species == "OHX(31)": 
#                 xmax = -1.0
#                 xmin = -4.5
                
#             elif species == "OX(25)": 
#                 xmax = -3.5 
#                 xmin = -6.5
                
#             # plot dft, lsr, ocp as the columns
#             axes = []
#             for col, src in enumerate(self.project_names):
#                 initval = self.spec_start_h298s[src][species]
#                 mapval = self.params_map_scaled[src][species]
#                 muapval = self.params_muap_scaled[src][species] 

#                 # our prior uncertainties
#                 ub = self.spec_uncs_ub[src][species]
#                 lb = self.spec_uncs_lb[src][species]
                
#                 # the posterior uncertainties
#                 ubstdap = mapval + stdaps_dict[src][species]/ev2jkm
#                 lbstdap = mapval - stdaps_dict[src][species]/ev2jkm
                
#                 ax = plt.Subplot(fig, inner[col])
#                 axes.append(ax)
#                 show_lit_vals = True
                    
#                 if show_lit_vals: 
#                     if species != "OX(25)" and species != "OHX(31)":
#                         ax.axvline(x=initval, color='#00A5DF', label = f"initial")
#                         ax.axvline(x=mapval, color='r', label = f"map")
#                         ax.axvline(x=muapval, color='g', label = f"$\mu$ap")
                    
#                     # source: abild pedersen paper
#                     elif species == "OHX(31)":
#                         ax.axvline(x=-3.6, color='r', label = f"literature (211)")

#                     # source: abild pedersen paper
#                     elif species == "OX(25)":
#                         ax.axvline(x=-5.0, color='r', label = f"literature (211)")
#                 else: 
#                     ax.axvline(x=initval, color='#00A5DF', label = f"initial")
#                     ax.axvline(x=mapval, color='r', label = f"map")
#                     ax.axvline(x=muapval, color='g', label = f"$\mu$ap")
                
                
#                 # make the xmin/xmax the max and min of the parameter uncertainty_vals
#                 ax.set_xlim(left=xmin, right=xmax)
                
#                 # plot normal distribution 
#                 mu = initval
#                 sigma = stdpriors_dict[src][species]/ev2jkm
#                 # sigma = stdaps_dict[src][species]/ev2jkm
#                 x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
#                 ax.fill_between(x, stats.norm.pdf(x, mu, sigma), alpha=0.2, facecolor="black", color=None, label="prior distribution" )
                
#                 # plot distribution of the post burn in data: 
#                 x_data = self.burn_in_bes[src][:,row]/ev2jkm + self.spec_start_h298s[src][species]
#                 ax.hist(x=x_data, bins = 1000, density=True, alpha = 0.2, color= "m", label="posterior distribution",) 
                
                
#                 # create a secondary y axis 
#                 ax2 = ax.twinx()
#                 ax2.set_ylabel(src,rotation='horizontal')
#                 ax2.set_yticklabels([])
#                 ax2.set_yticks([])
#                 ax2.yaxis.set_label_coords(1.05, 0.7)
                    
#                 # only show y-axis label if it is second row
#                 if col == 0:
#                     ax.get_xaxis().set_visible(False)
#                     ax.get_yaxis().set_visible(False)
#                 elif col ==1: 
#                     ax.get_xaxis().set_visible(False)
#                     ax.set_ylabel(pretty_species)
#                     ax.get_yaxis().set_visible(True)
#                 elif col == 2:
#                     ax.get_xaxis().set_visible(True)
#                     ax.set_xlabel("Binding Energy (eV)")
#                     ax.get_yaxis().set_visible(False)
#                     ax.legend(bbox_to_anchor=(1.1, 3.5), loc="upper left")
                
#                 # hide the y axis
#                 ax.set_yticklabels([])
#                 ax.set_yticks([])
#                 fig.add_subplot(ax)

#             # save only the legend: 
#             if row == 0: 
#                 legend = ax.get_legend()
#                 bbox = mpl.transforms.Bbox([[497.08333333333337, -12.718749999999986], [732.7083333333334, 110.94791666666667]])
#                 fig.canvas.draw()
#                 bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                
#                 # add a small buffer to the legend bounding box, so we don't get the legend outline
#                 bbox.x0 = bbox.x0 + abs(bbox.x0)*0.01
#                 bbox.x1 = bbox.x1 - abs(bbox.x1)*0.01
#                 bbox.y0 = bbox.y0 + abs(bbox.y0)*0.05
#                 bbox.y1 = bbox.y1 - abs(bbox.y1)*0.05
#                 fig.savefig('legend.pdf', dpi="figure",bbox_inches=bbox)
                
#             if species == "OX(25)": 
#                 legend = ax.get_legend()
#                 bbox = mpl.transforms.Bbox([[497.08333333333337, -12.718749999999986], [732.7083333333334, 110.94791666666667]])
#                 fig.canvas.draw()
#                 bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                
#                 # add a small buffer to the legend bounding box, so we don't get the legend outline
#                 bbox.x0 = bbox.x0 + abs(bbox.x0)*0.01
#                 bbox.x1 = bbox.x1 - abs(bbox.x1)*0.01
#                 bbox.y0 = bbox.y0 + abs(bbox.y0)*0.05
#                 bbox.y1 = bbox.y1 - abs(bbox.y1)*0.05
                
#                 if plot_path is not None: 
#                     save_path_legend = os.path.join(plot_path, f"{species}_legend.pdf")
#                     save_path_legend_lit = os.path.join(plot_path, f"{species}_legend_lit.pdf")
#                     fig.savefig(save_path_legend, dpi="figure",bbox_inches=bbox)
#                     fig.savefig(save_path_legend_lit, dpi="figure",bbox_inches=bbox)
            
#             if plot_path is not None: 
#                 save_path = os.path.join(plot_path, f"{species}_{plot_name}.pdf")
#                 fig.savefig(save_path, bbox_inches="tight")

#             plts.append((fig,axes,inner))
        
#         return plts

    


#     def make_convergence_plots(self, plot_path=None, plot_name=None): 
#         """
#         make geweke and autocorrellation plots for dataset
#         """
#         plot_paths = {}
#         for src in self.project_names:
#             # if not os.path.exists(os.path.join(project_path, "pickles/mcmc_post_burn_in_statistics.pkl")):
#             pe_object = PEUQSE.load_PE_object(os.path.join(self.project_paths[src], "pe_object.dill"))
#             param_names = pe_object.UserInput.model['parameterNamesAndMathTypeExpressionsDict']
#             sample_shape = pe_object.post_burn_in_samples.shape
#             chain_samples = pe_object.post_burn_in_samples.reshape(sample_shape[0], 1, sample_shape[1])
#             print(chain_samples.shape)
#             PEUQSE.InverseProblem.calculateAndPlotConvergenceDiagnostics(
#                 chain_samples,
#                 param_names,
#                 pe_object.UserInput.scatter_matrix_plots_settings, 
#                 pe_object.UserInput.directories['graphs'],
#                 showFigure = False,
#                 symbol=src,
#             )
#             plot_paths[src] = pe_object.UserInput.directories['graphs']
         
#         labels = [["a)", "b)"], ["c)", "d)"], ["e)", "f)"]]
#         fig, ax = plt.subplots(3,2)
#         fig.set_size_inches(3.5, 4)
#         # this determines the x-y position of the letters labelling each subplot
#         trans = mtransforms.ScaledTranslation(15/72, -15/72, fig.dpi_scale_trans)

#         for idx, (src, path) in enumerate(plot_paths.items()):
#             img_autocorr = plt.imread(os.path.join(path, f"AutoCorrelationPlot_Combined_Parameters.png"))
#             img_geweke = plt.imread(os.path.join(path, f"GewekeDiagnostic_Combined_Parameters.png"))
#             # print(src)
#             # display(img_autocorr)
#             # display(img_geweke)
#             ax[idx, 0].imshow(np.array(img_autocorr))
#             ax[idx, 1].imshow(np.array(img_geweke))

#             # hide the x and y axis and bounding box
#             ax[idx, 0].set_axis_off()
#             ax[idx, 0].set_aspect('equal')
#             ax[idx, 1].set_axis_off()
#             ax[idx, 1].set_aspect('equal')

#             ax[idx, 0].ylim = (-0.01, 1)

#             ax[idx, 0].text(0.0, 1.0, labels[idx][0], transform=ax[idx,0].transAxes + trans,
#                     fontsize=8, verticalalignment='top', #fontfamily='calibri',
#                     # bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0)
#             )
            
#             ax[idx, 1].text(0.0, 1.0, labels[idx][1], transform=ax[idx,1].transAxes + trans,
#                 fontsize=8, verticalalignment='top', #fontfamily='calibri',
#                 # bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0)
#             )

#         # make the subplots bigger and less whitespace between
#         fig.subplots_adjust(hspace=0.00, wspace=0.00, bottom=0.1, top=0.9, left=0.125, right=0.5)
#         fig.tight_layout()
        
#         if not plot_name: 
#             plot_name = "convergence_plots"
#         fig.savefig(os.path.join(plot_path, f"{plot_name}.pdf"), dpi=500,bbox_inches='tight')

#     def make_posterior_plots(self, plot_path=None, plot_name=None):
       
#         # pick which lit values we want to plot with each species
#         plotting_options = {
#             'HX(21)': '111',
#             'CO2X(22)': '111',
#             'COX(23)': '111',
#             'CH4X(24)': '111',
#             'OX(25)': '211',
#             'CH2X(26)': '111',
#             'CH3X(27)': '111',
#             'CHX(28)': '111',
#             'CX(29)': '111',
#             'H2X(30)': '111',
#             'OHX(31)': '211',
#             'H2OX(32)': '111',
#             'CHOX(33)': '111',
#             }

#         # import matplotlib.gridspec as gridspec
#         # import matplotlib.pyplot as plt

#         if not plot_path:
#             plot_path = os.path.join(repo_dir, "plotting", "plots")
#         else: 
#             plot_path = plot_path

#         # figure size
#         plt.rcParams["figure.figsize"] = [10,10]

#         # set font size
#         plt.rcParams.update({'font.size': 16})

#         # list of fig ax objects to return after, for use in jupyter
#         plts = []

#         # grab posterior stdev values
#         stdaps_dict = {}
#         stdpriors_dict = {}
#         for src in self.project_names:
#             stdaps_dict[src] = dict(zip(self.param_dicts[src]["label_list"], self.parameters_stdaps[src]))
#             stdpriors_dict[src] = dict(zip(self.param_dicts[src]["label_list"], self.parameters_stdpriors[src]))
            
#         outer = gridspec.GridSpec(len(self.param_dicts["DFT"]["label_list"]), 3, wspace=0.2, hspace=0.8)

#         fig = plt.figure(figsize=(20, 50))


#         self.species_list_pretty = [
#                     'H*',
#                     'CO$_2$*',
#                     'CO*',
#                     'CH$_4$*',
#                     'O* (211)',
#                     'CH$_2$*',
#                     'CH$_3$*',
#                     'CH*',
#                     'C*',
#                     'H$_2$*',
#                     'OH* (211)',
#                     'H$_2$O*',
#                     'CHO*'
#                 ]
#         self.species_dict = dict(zip(list(self.species_dict.keys()), self.species_list_pretty))
                                        
#         for row, (species, pretty_species) in enumerate(self.species_dict.items()):
#             # fig = plt.figure(figsize=(5, 3))
#             # inner = gridspec.GridSpec(3, 1, wspace=0.1, hspace=0.1)      
#             inner = gridspec.GridSpecFromSubplotSpec(3, 1,
#                     subplot_spec=outer[row], wspace=0.1, hspace=0.1)   

#             # subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)           
            
#             xmaxs = []
#             xmins = []

#             # define the max and min for the xaxis 
#             for src in self.project_names:
#                 xmaxs.append(self.spec_uncs_ub[src][species])
#                 xmins.append(self.spec_uncs_lb[src][species])
                
#                 # also append the max values from map, muap
#                 xmaxs.append(self.params_map_scaled[src][species])
#                 xmins.append(self.params_map_scaled[src][species])
            
#                 xmaxs.append(self.params_muap_scaled[src][species])
#                 xmins.append(self.params_muap_scaled[src][species])
                
#                 # also need to look at the new uncertainties for 
#                 # upper and lower bounds for the plot
#                 xmaxs.append(self.params_map_scaled[src][species] + 2*stdpriors_dict[src][species]/ev2jkm)
#                 xmins.append(self.params_map_scaled[src][species] - 2*stdpriors_dict[src][species]/ev2jkm)
            
#             # make the buffer 2*sigma so we can see the distribution better
#             # so far ocp has the largest uncertainty, So just hardcoding that
#             buffer = 2*stdaps_dict["ML"][species]/ev2jkm*2 
            
#             # # hardcoding max/min for oh and o
#             # if pretty_species != "OHX (211)" and pretty_species !="OX (211)": 
#             xmax = max(xmaxs) + buffer
#             xmin = min(xmins) - buffer
                
#             # elif pretty_species == "OHX (211)": 
#             #     xmax = -1.0
#             #     xmin = -4.5
                
#             # elif pretty_species == "OX (211)": 
#             #     xmax = -3.5 
#             #     xmin = -6.5
                
#             # plot dft, lsr, ocp as the columns
#             axes = []
#             for col, src in enumerate(self.project_names):
#                 initval = self.spec_start_h298s[src][species]
#                 mapval = self.params_map_scaled[src][species]
#                 muapval = self.params_muap_scaled[src][species] 

#                 # our prior uncertainties
#                 ub = self.spec_uncs_ub[src][species]
#                 lb = self.spec_uncs_lb[src][species]
                
#                 # the posterior uncertainties
#                 ubstdap = mapval + stdaps_dict[src][species]/ev2jkm
#                 lbstdap = mapval - stdaps_dict[src][species]/ev2jkm
                
#                 ax = plt.Subplot(fig, inner[col])
#                 axes.append(ax)
#                 show_lit_vals = False
                    
#                 if show_lit_vals: 

#                     # change evaluation to if lit values exist or not
#                     lit_be_111 = lit_vals_111[species]
#                     lit_be_211 = lit_vals_211[species]
#                     if lit_be_111 and plotting_options[species] == "111":
#                         for lit_src in lit_be_111.keys():
#                             ax.axvline(x=lit_be_111[lit_src], color='b', label = f"Literature value 111/211")
#                     if lit_be_211 and plotting_options[species] == "211":
#                         for coloridx,lit_src in enumerate(lit_be_211.keys()):
#                             if lit_src == "zhang2011":
#                                 ax.axvline(x=lit_be_211[lit_src], color='b', label = f"Literature value 111/211")

#                     # if pretty_species != "OHX (211)" and pretty_species != "OX (211)":
#                         # ax.axvline(x=initval, color='#00A5DF', label = f"initial")
#                         # ax.axvline(x=mapval, color='r', label = f"map")
#                         # ax.axvline(x=muapval, color='g', label = f"$\mu$ap")
#                     # draw a vertical line for the map value
                            
#                     ax.axvline(x=mapval, color='r', label = f"map")


#                     # # source: abild pedersen paper
#                     # elif pretty_species == "OHX (211)":
#                     #     ax.axvline(x=-3.6, color='r', label = f"literature (211)")

#                     # # source: abild pedersen paper
#                     # elif pretty_species == "OX (211)":
#                     #     ax.axvline(x=-5.0, color='r', label = f"literature (211)")
                    
#                 else: 
#                     # ax.axvline(x=initval, color='#00A5DF', label = f"initial")
#                     ax.axvline(x=mapval, color='r', label = f"map")
#                     # ax.axvline(x=muapval, color='g', label = f"$\mu$ap")
                
                
#                 # make the xmin/xmax the max and min of the parameter uncertainty_vals
#                 ax.set_xlim(left=xmin, right=xmax)
                
#                 # plot normal distribution 
#                 mu = initval
#                 # sigma = stdaps_dict[src][species]/ev2jkm
#                 sigma = stdpriors_dict[src][species]/ev2jkm
#                 x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
#                 ax.fill_between(x, stats.norm.pdf(x, mu, sigma), alpha=0.2, facecolor="black", color=None, label="prior distribution" )
                
#                 # plot distribution of the post burn in data: 
#                 x_data = self.burn_in_bes[src][:,row]/ev2jkm + self.spec_start_h298s[src][species]
#                 ax.hist(x=x_data, bins = 1000, density=True, alpha = 0.2, color= "m", label="posterior distribution",) 
                
                
#                 # create a secondary y axis 
#                 ax2 = ax.twinx()
#                 ax2.set_ylabel(src,rotation='horizontal')
#                 ax2.set_yticklabels([])
#                 ax2.set_yticks([])
#                 ax2.yaxis.set_label_coords(1.05, 0.7)
                    
#                 # only show y-axis label if it is second row
#                 if col == 0:
#                     ax.get_xaxis().set_visible(False)
#                     ax.get_yaxis().set_visible(False)
#                 elif col ==1: 
#                     ax.get_xaxis().set_visible(False)
#                     ax.set_ylabel(pretty_species)
#                     ax.get_yaxis().set_visible(True)
#                 elif col == 2:
#                     ax.get_xaxis().set_visible(True)
#                     ax.set_xlabel("Binding Energy (eV)")
#                     ax.get_yaxis().set_visible(False)
#                     # ax.legend(bbox_to_anchor=(1.1, 3.5), loc="upper left")
                
#                 # hide the y axis
#                 ax.set_yticklabels([])
#                 ax.set_yticks([])

#                 # show the legend if it's the last plot 
#                 if row == len(self.species_dict.keys())-1 and col == 2:
#                     ax.legend(bbox_to_anchor=(2.0, 1.5), loc="center right")

#                 fig.add_subplot(ax)
        
#         if not plot_name: 
#             plot_name = "all_species_prob_density"
            
#         save_path = os.path.join(plot_path, f"{plot_name}.pdf")
#         fig.savefig(save_path, bbox_inches="tight")
        

#     def make_combined_output_plots(self, plot_path=None, plot_name=None, no_exp_dat_path=None, alpha = 0.1): 
#         # plot most probable and best fit set
#         plt.rcParams["figure.figsize"] = [20,25]
#         fig, ax = plt.subplots(len(self.value_list)-1,3)
#         # fig, ax = plt.subplots(1,3, sharey=True)

#         y_axis_labels = {
#             'ch4_profile': 'CH$_4$',
#             'co2_profile': 'CO$_2$',
#             'co_profile': 'CO',
#             'h2_profile': 'H$_2$',
#             'h2o_profile': 'H$_2$O',
#             'o2_profile': 'O$_2$'
#             }

#         plts = []
#         for axis, op in enumerate(self.value_list):
#             if self.label_list[axis] != "temp_profile":
#                 # fig, ax = plt.subplots(1,3, sharey=True)

#                 # loop over lsr, dft, ocp results
#                 for idx, src in enumerate(self.project_names):
#                     ax[axis,idx].scatter(self.x_vals, self.output_initials[src][axis], label=f"Initial", color="g")
#                     ax[axis,idx].scatter(self.x_vals, self.output_maps[src][axis], label=f"MAP", color="orange")
#                     # ax[axis,idx].scatter(self.x_vals, self.output_mu_aps[src][axis], label=f"{src} - mu_ap", color="orange")
#                     # ax[axis,idx].scatter(self.x_vals, self.output_std_aps[src][axis], label=f"{src} - std_ap")
#                     ax[axis,idx].plot(self.x_vals, op, label="Experiment", linestyle="-", linewidth=1.5, marker="o")
#                     # ax[axis,idx].errorbar(self.x_vals, op, yerr=self.uncerts[self.label_list[axis]], capsize=2.5, capthick=1)
#                     # ax[axis].errorbar(self.x_vals, op, yerr=flat_self.uncerts[self.label_list[axis]], capsize=2.5, capthick=1)
#                     # ax[axis].errorbar(self.x_vals, op, yerr=self.uncerts[self.label_list[axis]], capsize=2.5, capthick=1)
#                     # ax[axis,idx].errorbar(self.x_vals, self.output_mu_aps[src][axis], yerr=self.burn_in_outputs_std[axis], capsize=2.5, capthick=1)

#                     # trim lower error bars if they go below 0
#                     exp_upper_error = np.array(op) + np.array(self.uncerts[self.label_list[axis]])*2
#                     exp_lower_error = np.array(op) - np.array(self.uncerts[self.label_list[axis]])*2
#                     exp_lower_error[exp_lower_error < 0] = 0

#                     sim_upper_error = self.output_maps[src][axis] + self.burn_in_outputs_std[src][axis]*2
#                     sim_lower_error = self.output_maps[src][axis] - self.burn_in_outputs_std[src][axis]*2
#                     sim_lower_error[sim_lower_error < 0] = 0


#                     ax[axis,idx].fill_between(
#                         self.x_vals, 
#                         # self.output_mu_aps[src][axis] - self.burn_in_outputs_std[axis],
#                         # self.output_mu_aps[src][axis] + self.burn_in_outputs_std[axis],
#                         sim_lower_error,
#                         sim_upper_error,
#                         color="sandybrown",
#                         alpha= 0.3,
#                         label="2$\sigma$ MAP",
#                         )
#                     ax[axis,idx].fill_between(
#                         self.x_vals, 
#                         exp_lower_error,
#                         exp_upper_error,
#                         color="powderblue",
#                         alpha= 0.3,
#                         label="2$\sigma$ Experiment"
#                         )

#                     # if we specify it, include uncertainty for initial data when we exclude experimental data uncertainty
#                     if no_exp_dat_path: 
#                         with open(no_exp_dat_path, "rb") as f:
#                             initial_point_err = pickle.load(f)
                            
#                         init_upper_error = self.output_initials[src][axis] + initial_point_err[src][axis]*2
#                         init_lower_error = self.output_initials[src][axis] - initial_point_err[src][axis]*2
#                         init_lower_error[init_lower_error < 0] = 0
#                         ax[axis,idx].fill_between(
#                             self.x_vals, 
#                             init_upper_error,
#                             init_lower_error,
#                             color="g",
#                             alpha=alpha,
#                             label="2$\sigma$ Initial"
#                             )

#                     if axis == 2 and idx == 2: 
#                         ax[axis,idx].legend(bbox_to_anchor=(1.2,0.9))

#                     if axis == 0:
#                         ax[axis,idx].set_title(f"{src} ")
#                     # print(len(self.value_list)-1, axis)
#                     if axis == len(self.value_list)-2:
#                         ax[axis,idx].get_xaxis().set_visible(True)
#                         ax[axis,idx].set_xlabel("Reactor coordinate (mm)")
#                     else: 
#                          ax[axis,idx].get_xaxis().set_visible(False)

#                     if idx == 0:
#                         ax[axis,idx].set_ylabel(f"{y_axis_labels[self.label_list[axis]]} outlet flow (mol/min)")
#                     else: 
#                         ax[axis,idx].get_yaxis().set_visible(False)


#                     if idx ==2: 
#                         # ax[axis,idx].legend(bbox_to_anchor=(1.5,0.9))
#                         buff = 0.15
#                         # at end, readjust the y limits to make sure they are all identical
#                         ymax = max([ax[axis,0].get_ylim()[1], ax[axis,1].get_ylim()[1], ax[axis,2].get_ylim()[1]])
#                         ymin = min([ax[axis,0].get_ylim()[0], ax[axis,1].get_ylim()[0], ax[axis,2].get_ylim()[0]])
#                         ax[axis,0].set_ylim(ymin, ymax+ymax*buff)
#                         ax[axis,1].set_ylim(ymin, ymax+ymax*buff)
#                         ax[axis,2].set_ylim(ymin, ymax+ymax*buff)

#                     fig.canvas.draw()

#         if not plot_name: 
#             plot_name = "all_output_plots"
#         if plot_path is not None: 
#             fig.subplots_adjust(wspace=0.1, hspace=0.1)
#             save_path = os.path.join(plot_path, f"{plot_name}.pdf")
#             fig.savefig(save_path,bbox_inches="tight")
            
            
            
#     def make_act_plots(self, plot_path=None, plot_name=None):
        
#         plt.rcParams["figure.figsize"] = [15,5]
#         fig, ax = plt.subplots(1,3)

#         for col, src in enumerate(self.project_names): 
#             pe_object = PEUQSE.load_PE_object(os.path.join(self.project_paths[src], "pe_object.dill"))
#             param_names = pe_object.UserInput.model['parameterNamesAndMathTypeExpressionsDict']
#             sample_shape = pe_object.post_burn_in_samples.shape
#             chain_samples = pe_object.post_burn_in_samples.reshape(sample_shape[0], 1, sample_shape[1])
#             plot_outputs = PEUQSE.InverseProblem.calculateAndPlotConvergenceDiagnostics(
#                 chain_samples,
#                 param_names,
#                 pe_object.UserInput.scatter_matrix_plots_settings, 
#                 pe_object.UserInput.directories['graphs'],
#                 showFigure = False,
#                 symbol=src,
#             )
#             # some of the initial values are nan or <0, so skip those
#             start_id = 4
#             ax[col].loglog(plot_outputs[0][start_id:], plot_outputs[0][start_id:]/50,"--k", label=r"$\tau = N/50$")
            
#             for key, val in plot_outputs[2].items():
#                 # only plot values > a threshold value (lsr plot has an initial value that is really small 
#                 # and throws off plot
#                 plot_ids = np.where(val>1e-10)
#                 ax[col].loglog(plot_outputs[0][plot_ids], val[plot_ids], marker="", linestyle="-", label=self.species_dict[key])
#             if col == 0: 
#                 ax[col].set(ylabel=r"$\tau$ estimates")
#             if col > 0: 
#                 ax[col].get_yaxis().set_visible(False)
#             if col == 2: 
#                 ax[col].legend(bbox_to_anchor=(1.1, 1.0))

#             ax[col].set(xlabel="number of samples, $N$")

#             ax[col].set_title(f"{src} ")
        
#         if not plot_name: 
#             plot_name = "act_plots"
#         fig.savefig(os.path.join(plot_path, f"{plot_name}.pdf"),bbox_inches = 'tight')\
        
        
#     def plot_one_be(self, species, src, fig=None, ax = None, plot_path=None, plot_name=None):
            
#         """ 
#         plot the posterior distribution for the molar flow rate data
#         """
        
#         if not plot_name: 
#             plot_name = "figure"
#         # load the 
#         # figure size
#         plt.rcParams["figure.figsize"] = [15,3]

#         # set font size
#         plt.rcParams.update({'font.size': 12})

#         # pick which lit values we want to plot with each species
#         plotting_options = {
#             'HX(21)': '111',
#             'CO2X(22)': '111',
#             'COX(23)': '111',
#             'CH4X(24)': '111',
#             'OX(25)': '211',
#             'CH2X(26)': '111',
#             'CH3X(27)': '111',
#             'CHX(28)': '111',
#             'CX(29)': '111',
#             'H2X(30)': '111',
#             'OHX(31)': '211',
#             'H2OX(32)': '111',
#             'CHOX(33)': '111',
#             }
        
#         # grab posterior stdev values
#         stdaps_dict = {}
#         stdpriors_dict = {}

#         stdaps_dict[src] = dict(zip(self.param_dicts[src]["label_list"], self.parameters_stdaps[src]))
#         stdpriors_dict[src] = dict(zip(self.param_dicts[src]["label_list"], self.parameters_stdpriors[src]))
            
#         row = list(self.species_dict.keys()).index(species)
#         pretty_species = self.species_dict[species]

#         self.species_list_pretty = [
#                     'H*',
#                     'CO$_2$*',
#                     'CO*',
#                     'CH$_4$*',
#                     'O* (211)',
#                     'CH$_2$*',
#                     'CH$_3$*',
#                     'CH*',
#                     'C*',
#                     'H$_2$*',
#                     'OH* (211)',
#                     'H$_2$O*',
#                     'CHO*'
#                 ]
#         self.species_dict = dict(zip(list(self.species_dict.keys()), self.species_list_pretty))

#         xmaxs = []
#         xmins = []

#         # define the max and min for the xaxis 
#         for src in self.project_names:
#             xmaxs.append(self.spec_uncs_ub[src][species])
#             xmins.append(self.spec_uncs_lb[src][species])

#             # also append the max values from map, muap
#             xmaxs.append(self.params_map_scaled[src][species])
#             xmins.append(self.params_map_scaled[src][species])

#             xmaxs.append(self.params_muap_scaled[src][species])
#             xmins.append(self.params_muap_scaled[src][species])

#             # also need to look at the new uncertainties for 
#             # upper and lower bounds for the plot
#             xmaxs.append(self.params_map_scaled[src][species] + 2*stdpriors_dict[src][species]/ev2jkm)
#             xmins.append(self.params_map_scaled[src][species] - 2*stdpriors_dict[src][species]/ev2jkm)

#         # make the buffer 2*sigma so we can see the distribution better
#         # so far ocp has the largest uncertainty, So just hardcoding that
#         buffer = 2*stdaps_dict["ML"][species]/ev2jkm*2 

#         # # hardcoding max/min for oh and o
#         # if pretty_species != "OHX (211)" and pretty_species !="OX (211)": 
#         xmax = max(xmaxs) + buffer
#         xmin = min(xmins) - buffer

#         # plot dft, lsr, ocp as the columns
#         initval = self.spec_start_h298s[src][species]
#         mapval = self.params_map_scaled[src][species]
#         muapval = self.params_muap_scaled[src][species] 

#         # our prior uncertainties
#         ub = self.spec_uncs_ub[src][species]
#         lb = self.spec_uncs_lb[src][species]

#         # the posterior uncertainties
#         ubstdap = mapval + stdaps_dict[src][species]/ev2jkm
#         lbstdap = mapval - stdaps_dict[src][species]/ev2jkm

#         # ax = plt.Subplot(fig, inner[col])
#         axes.append(ax)
#         show_lit_vals = False

#         ax.axvline(x=mapval, color='r', label = f"map")

#         # make the xmin/xmax the max and min of the parameter uncertainty_vals
#         ax.set_xlim(left=xmin, right=xmax)

#         # plot normal distribution 
#         mu = initval
#         sigma = stdpriors_dict[src][species]/ev2jkm
#         x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
#         ax.fill_between(x, stats.norm.pdf(x, mu, sigma), alpha=0.2, facecolor="black", color=None, label="prior distribution" )

#         # plot distribution of the post burn in data: 
#         x_data = self.burn_in_bes[src][:,row]/ev2jkm + self.spec_start_h298s[src][species]
#         ax.hist(x=x_data, bins = 1000, density=True, alpha = 0.2, color= "m", label="posterior distribution",) 

#         # create a secondary y axis 
#         ax2 = ax.twinx()
#         ax2.set_ylabel(src,rotation='horizontal')
#         ax2.set_yticklabels([])
#         ax2.set_yticks([])
#         ax2.yaxis.set_label_coords(1.05, 0.7)

#         # only show y-axis label if it is second row
#         if col == 0:
#             ax.get_xaxis().set_visible(False)
#             ax.get_yaxis().set_visible(False)
#         elif col ==1: 
#             ax.get_xaxis().set_visible(False)
#             ax.set_ylabel(pretty_species)
#             ax.get_yaxis().set_visible(True)
#         elif col == 2:
#             ax.get_xaxis().set_visible(True)
#             ax.set_xlabel("Binding Energy (eV)")
#             ax.get_yaxis().set_visible(False)
#             # ax.legend(bbox_to_anchor=(1.1, 3.5), loc="upper left")

#         # hide the y axis
#         ax.set_yticklabels([])
#         ax.set_yticks([])

#         # show the legend if it's the last plot 
#         if row == len(self.species_dict.keys())-1 and col == 2:
#             ax.legend(bbox_to_anchor=(2.0, 1.5), loc="center right")

#         fig.add_subplot(ax)

#         if not plot_name: 
#             plot_name = f"{species}_species_prob_density"

#         save_path = os.path.join(plot_path, f"{plot_name}.pdf")
#         fig.savefig(save_path, bbox_inches="tight")




       