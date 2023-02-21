###############################################################################
# preprocess.py
# preprocessing script for interpreting yaml input from experimental data. 
# Chris Blais
# northeastern university
# 2022
###############################################################################

import yaml
import sys
import itertools
import logging


'''
arguments to script: 
# python preprocess.py <yaml_file> <expt_type>
# yaml_file = yaml file containing experimental conditions
# expt_type: either "sbr" or "batch" for whether it is a single crystal batch 
# reactor or a spinning basket reactor.
'''
yaml_expt_file = sys.argv[1] #'all_experiments_temp_study.yaml'
expt_type = sys.argv[2] #"sbr"

def zip_dict(**kwargs):
    '''
    a tool for iterating through a dictionary where the keys are the variable 
    (temp, pressure, etc) and the values are a list of possible values. 
    copied from stack overflow: 
    https://tinyurl.com/4erzdxeb
    '''
    keys = kwargs.keys()
    vals = kwargs.values()
    for instance in itertools.zip_longest(*vals):
        yield dict(zip(keys, instance))

with open(yaml_expt_file) as f:
    # use safe_load instead load
    expt_map = yaml.safe_load(f)
    
expt_list = []
for expt in expt_map.keys():
    if all(type_e == expt_type for type_e in expt_map[expt]['experiment_type']):
        expt_list.append(expt)

# list that retains expt variable names, 
# with a single value, and saves as a dictionary
expt_dict_list =[]
for name in expt_list:
    if not expt_dict_list:
        expt_dict_list = list(zip_dict(**expt_map[name]))
    else: 
        expt_dict_list.extend(list(zip_dict(**expt_map[name])))

yaml_reorg = "all_experiments_reorg_" + expt_type + ".yaml"

with open(yaml_reorg, 'w') as f:
    # use safe_load instead load
    doc = yaml.safe_dump(expt_dict_list, f)

# make another yaml with just the opt expts 
with open("./all_experiments_reorg_sbr.yaml", "r") as f:
    expt_all = yaml.load(f, Loader=yaml.FullLoader)

opt_expt_file = "./experiments_reorg_onlyopt.yaml"
expt_copy = []
for expt in expt_all:
    if "use_for_opt" in expt.keys(): 
        expt_copy.append(expt)

with open(opt_expt_file, "w") as f:
    yaml.dump(expt_copy, f)