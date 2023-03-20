# runtime utilities
# these functions are strictly for running peuqse. since peuqse is incompatible with rmg env, 
# these must only use peuqse modules. it's a mess but until rmg env is less static this is 
# how it has to be done. 

import sys 
import os
import yaml

def get_all_param_lists(results_path, kinetic=True, thermo=True):
    """
    returns lists for all the parameters used in the peuqse model
    kinetic - get all kinetic parameters (A, E0, alpha from rmg rules)
    thermo - get all thermo parameters (CHON and vdw be)
    if kinetic and thermo are true, append thermo parameters to kinetic parameters

    returns:
    value list - list of all parameter values
    label list - list of all parameter labels
    unc_list - list of all parameter uncertainties
    upper list - list of all parameter upper bounds
    lower list - list of all parameter lower bounds
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
    with open(os.path.join(results_path, "be_values.yaml"), 'r') as f:
        thermo_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "be_unc.yaml"), 'r') as f:
        thermo_unc_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "be_lb.yaml"), 'r') as f:
        thermo_lb_config = yaml.safe_load(f)
    with open(os.path.join(results_path, "be_ub.yaml"), 'r') as f:
        thermo_ub_config = yaml.safe_load(f)

    value_list = []
    label_list = []
    unc_list = []
    upper_list = []
    lower_list = []

    if kinetic and thermo: 
        # kinetic parameters
        for label, value in rule_config.items():
            value_list.extend([value["A"], value["E0"], value["alpha"]])
            label_list.extend([label + " A", label + " E0", label + " alpha"])
            unc_list.extend([rule_unc_config[label]["A"], rule_unc_config[label]["E0"], rule_unc_config[label]["alpha"]])
            upper_list.extend([rule_ub_config[label]["A"], rule_ub_config[label]["E0"], rule_ub_config[label]["alpha"]])
            lower_list.extend([rule_lb_config[label]["A"], rule_lb_config[label]["E0"], rule_lb_config[label]["alpha"]])

        # thermo parameters
        for label, value in thermo_config.items():
            value_list.append(value)
            label_list.append(label)
            unc_list.append(thermo_unc_config[label])
            upper_list.append(thermo_ub_config[label])
            lower_list.append(thermo_lb_config[label])
    # return as a dictionary so we don't confuse what's what
    peuqse_params = {
        "value_list": value_list,
        "label_list": label_list,
        "unc_list": unc_list,
        "upper_list": upper_list,
        "lower_list": lower_list
    }
    return peuqse_params
