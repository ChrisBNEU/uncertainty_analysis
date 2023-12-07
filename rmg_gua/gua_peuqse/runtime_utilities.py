# runtime utilities
# these functions are strictly for running peuqse. since peuqse is incompatible with rmg env, 
# these must only use peuqse modules. it's a mess but until rmg env is less static this is 
# how it has to be done. 

import sys 
import os
import yaml

def get_all_param_lists(results_path, kinetic=True, thermo=True, reduce_space=None):
    """
    returns lists for all the parameters used in the peuqse model
    kinetic - get all kinetic parameters (A, E0, alpha from rmg rules)
    thermo - get all thermo parameters (CHON and vdw be)
    if kinetic and thermo are true, append thermo parameters to kinetic parameters

    reduce_space - (optional) remove parameters from peuqse run, by specifying
    the indices of the parameters that we explicitly want to keep, or by specifying 
    a string we want to match. 

    returns:
    value list - list of all parameter values
    label list - list of all parameter labels
    unc_list - list of all parameter uncertainties
    upper list - list of all parameter upper bounds
    lower list - list of all parameter lower bounds
    reduce_list - (optional) list of the indices for the parameters you want to run
    """
    # open all of the config files
    with open(os.path.join(results_path, "rule_config.yaml"), 'r') as f:
        rule_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "rule_unc_config.yaml"), 'r') as f:
        rule_unc_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "rule_lb_config.yaml"), 'r') as f:
        rule_lb_config = yaml.safe_load(f) 
    with open(os.path.join(results_path, "rule_ub_config.yaml"), 'r') as f:
        rule_ub_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "rule_guess_config.yaml"), 'r') as f:
        rule_guess_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "be_values.yaml"), 'r') as f:
        thermo_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "be_unc.yaml"), 'r') as f:
        thermo_unc_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "be_lb.yaml"), 'r') as f:
        thermo_lb_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "be_ub.yaml"), 'r') as f:
        thermo_ub_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "be_guess.yaml"), 'r') as f:
        thermo_guess_config = yaml.safe_load(f)

    value_list = []
    label_list = []
    unc_list = []
    upper_list = []
    lower_list = []
    guess_list = []

    if kinetic and thermo: 
        # kinetic parameters
        for label, value in rule_config.items():
            value_list.extend([value["A"], value["E0"], value["alpha"]])
            label_list.extend([label + " A", label + " E0", label + " alpha"])
            unc_list.extend([rule_unc_config[label]["A"], rule_unc_config[label]["E0"], rule_unc_config[label]["alpha"]])
            upper_list.extend([rule_ub_config[label]["A"], rule_ub_config[label]["E0"], rule_ub_config[label]["alpha"]])
            lower_list.extend([rule_lb_config[label]["A"], rule_lb_config[label]["E0"], rule_lb_config[label]["alpha"]])
            guess_list.extend([rule_guess_config[label]["A"], rule_guess_config[label]["E0"], rule_guess_config[label]["alpha"]])
                               
        # thermo parameters
        for label, value in thermo_config.items():
            value_list.append(value)
            label_list.append(label)
            unc_list.append(thermo_unc_config[label])
            upper_list.append(thermo_ub_config[label])
            lower_list.append(thermo_lb_config[label])
            guess_list.append(thermo_guess_config[label])
    
    # only do thermo parameter perturbation. more for testing 
    elif not kinetic and thermo:
        for label, value in thermo_config.items():
            value_list.append(value)
            label_list.append(label)
            unc_list.append(thermo_unc_config[label])
            upper_list.append(thermo_ub_config[label])
            lower_list.append(thermo_lb_config[label])
            guess_list.append(thermo_guess_config[label])
        
    # return as a dictionary so we don't confuse what's what
    peuqse_params = {
        "value_list": value_list,
        "label_list": label_list,
        "unc_list": unc_list,
        "upper_list": upper_list,
        "lower_list": lower_list,
        "guess_list": guess_list, 
    }

    reduce_dict = {}
    if reduce_space is not None: 
        # remove parameters from peuqse run, by specifying the indices of the 
        # parameters that we explicitly want to keep.
        print("reducing space")
        if isinstance(reduce_space, list):
            # simply supply the exact indices you want to exclude
            if all(isinstance(entry, int) for entry in reduce_space):
                reduce_list = reduce_space

                for idx, label in enumerate(label_list): 
                    if idx not in reduce_space: 
                        reduce_dict[label] = idx
                    else:
                        reduce_dict[label] = "removed"
            
            # if it is a list of strings, remove the strings that end with our parameter
            elif all(isinstance(entry, str) for entry in reduce_space): 
                print("reducing list of strings")
                reduce_list = []
                for label in label_list: 
                    if all(not label.endswith(entry) for entry in reduce_space): 
                        reduce_list.append(label_list.index(label))
                        reduce_dict[label] = label_list.index(label)
                    else:
                        reduce_dict[label] = "removed"
            else: 
                raise Exception("please specify either a list of strings or integers for reducing params")

        elif isinstance(reduce_space, str):
            # supply a string match for the parameters we want to exclude. 
            reduce_list = []

            for label in label_list: 
                if reduce_space not in label: 
                    reduce_list.append(label_list.index(label))
                    reduce_dict[label] = label_list.index(label)
                else: 
                    reduce_dict[label] = "removed"

        peuqse_params["reduce_list"] = reduce_list

        # save a yaml for the reduced list as well
        with open(os.path.join(results_path, "reduce_list.yaml"), 'w') as f:
            yaml.dump(reduce_dict, f, sort_keys=False)
    

    return peuqse_params


def make_exp_data_lists(results_path, use_count=False):
    """ 
    make a list of lists for the experimental data from graaf. 
    each sublist is an experiment. 
    for saving/loading yamls, we want to not sort yaml dictionaries. 
    input is expected to be a list of dictionaries, with each dictionary 
    having key/value pairs for temp, pressure, etc. 
    """

    expt_yaml_path = os.path.join(results_path, "expt_data.yaml")
    with open(expt_yaml_path, 'r') as f:
        expt_yaml = yaml.safe_load(f)

    expt_yaml_path = os.path.join(results_path, "expt_unc.yaml")
    with open(expt_yaml_path, 'r') as f:
        expt_unc_yaml = yaml.safe_load(f)
    
    # make the list of lists.
    # x_data is expected to have all experimental values in each sub list.
    # e.g [[T1, P1, V1, ...], [T2, P2, V2, ...], ... to nexpts] 
    # y_data is expected to have each sub list be on experimental value
    # e.g [[X1, X2, X3, ...], [Y1, Y2, Y3, ... ], ... to n_output parameters]
    if use_count:
        fcount = use_count
        
    
    exp_ct_input = []
    exp_list_in = []
    exp_list_in_alt = [[],[],[],[],[],[],[],]
    exp_list_y = [[],[],[],[],[],]
    exp_unc_list_y = [[],[],[],[],[],]
    count = 0


    exp_list_in = [
        expt_yaml["catalyst_area"],
        expt_yaml["pressure"],
        expt_yaml["species_CO"],
        expt_yaml["species_CO2"],
        expt_yaml["species_H2"],
        expt_yaml["temperature"],
        expt_yaml["volume_flowrate"]
        ]


    exp_list_y = [
        expt_yaml['species_out_CH3OH'],
        expt_yaml['species_out_CO'],
        expt_yaml['species_out_CO2'],
        expt_yaml['species_out_H2'],
        expt_yaml['species_out_H2O'],
    ]
    
    
    exp_unc_list_y = [
        expt_unc_yaml['species_out_CH3OH'],
        expt_unc_yaml['species_out_CO'],
        expt_unc_yaml['species_out_CO2'],
        expt_unc_yaml['species_out_H2'],
        expt_unc_yaml['species_out_H2O'],
    ]
    
    return exp_list_in, exp_list_y, exp_unc_list_y
        