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


def make_exp_data_lists(data_path, results_path):
    """ 
    make a list of lists for the experimental data from graaf. 
    each sublist is an experiment. 
    for saving/loading yamls, we want to not sort yaml dictionaries. 
    input is expected to be a list of dictionaries, with each dictionary 
    having key/value pairs for temp, pressure, etc. 
    """

    expt_yaml_path = os.path.join(data_path, "all_experiments_reorg_sbr.yaml")
    with open(expt_yaml_path, 'r') as f:
        expt_yaml = yaml.safe_load(f)

    # make the list of lists.
    # x_data is expected to have all experimental values in each sub list.
    # e.g [[T1, P1, V1, ...], [T2, P2, V2, ...], ... to nexpts] 
    # y_data is expected to have each sub list be on experimental value
    # e.g [[X1, X2, X3, ...], [Y1, Y2, Y3, ... ], ... to n_output parameters]

    exp_list_in = []
    expt_list_y = [[],[],[],[],[],]
    expt_unc_list_y = [[],[],[],[],[],]
    for exp in expt_yaml:
        exp_list_in.append(
            exp["catalyst_area"],
            exp["pressure"],
            exp["species"]["CO"],
            exp["species"]["CO2"],
            exp["species"]["H2"],
            exp["temperature"],
            exp["volume_flowrate"],
            )
        
        exp_list_y[0].append(exp['species_out']['CH3OH'])
        exp_list_y[1].append(exp['species_out']['CO'])
        exp_list_y[2].append(exp['species_out']['CO2'])
        exp_list_y[3].append(exp['species_out']['H2'])
        exp_list_y[4].append(exp['species_out']['H2O'])

        exp_unc_list_y[0].append(exp['species_out']['CH3OH']*0.01)
        exp_unc_list_y[1].append(exp['species_out']['CO']*0.01)
        exp_unc_list_y[2].append(exp['species_out']['CO2']*0.01)
        exp_unc_list_y[3].append(exp['species_out']['H2']*0.01)
        exp_unc_list_y[4].append(exp['species_out']['H2O']*0.01) 

    return exp_list_in, exp_list_y, exp_unc_list_y
        