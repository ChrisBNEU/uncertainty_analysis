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
    return peuqse_params


def make_exp_data_lists(data_path, results_path, use_count=False):
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
    if use_count:
        fcount = use_count
        
    
    exp_ct_input = []
    exp_list_in = []
    exp_list_in_alt = [[],[],[],[],[],[],[],]
    exp_list_y = [[],[],[],[],[],]
    exp_unc_list_y = [[],[],[],[],[],]
    count = 0
    for exp in expt_yaml:
        if "use_for_opt" in exp.keys() and exp["use_for_opt"]: 
            
            # append 
            exp_ct_input.append(exp)
            
            exp_list_in.append([
                exp["catalyst_area"],
                exp["pressure"],
                exp["species"]["CO"],
                exp["species"]["CO2"],
                exp["species"]["H2"],
                exp["temperature"],
                exp["volume_flowrate"]
            ])
            # trying alternate method because I am getting errors with the above
            exp_list_in_alt[0].append(exp["catalyst_area"])
            exp_list_in_alt[1].append(exp["pressure"])
            exp_list_in_alt[2].append(exp["species"]["CO"])
            exp_list_in_alt[3].append(exp["species"]["CO2"])
            exp_list_in_alt[4].append(exp["species"]["H2"])
            exp_list_in_alt[5].append(exp["temperature"])
            exp_list_in_alt[6].append(exp["volume_flowrate"])
                                      
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
            
            count+=1
            
        if use_count and count >=fcount:
            break
            
    
    exp_ct_file = os.path.join(results_path, "ct_expt_list.yaml")
    with open(exp_ct_file, "w") as f:
        yaml.safe_dump(exp_ct_input, f, sort_keys=False)
        

    return exp_list_in, exp_list_in_alt, exp_list_y, exp_unc_list_y
        