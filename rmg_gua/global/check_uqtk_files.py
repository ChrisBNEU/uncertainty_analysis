#!/usr/bin/env python
# coding: utf-8

# ## check the input values to our uncertainty run
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


# make pdf pages obj
pp = PdfPages('range_check.pdf')

# files for uqtk
pnamefile = "parameter_names_meoh.txt"
prangefile = "parameter_ranges_meoh.txt"
inpfile = "Input_meoh.txt"

param_names = []
param_ranges = []
inputs = []

input_dict = {}
range_dict = {}

with open(pnamefile, 'r') as f:
    for line in f.readlines():
        line_clean = line.replace("\n", "")
        param_names.append(line_clean)
        
with open(prangefile, 'r') as f:
    for line in f.readlines():
        values = line.split(" ")
        values = [float(val) for val in values]
        param_ranges.append(values)

with open(inpfile, 'r') as f:
    for line in f.readlines():
        values = line.split(" ")
        values = [float(val) for val in values]
        inputs.append(values) 
        
for num, val in enumerate(param_names):
    input_dict[val] = inputs[num]

for num, val in enumerate(param_names):
    inpt_list = [run[num] for run in inputs]
    inpt_list = np.array(inpt_list)
    input_dict[val] = inpt_list
    range_dict[val] = param_ranges[num]

# This just plots all of your ranges and values and ensures that 
# the perturbed values fall within the specified ranges (sanity check)
for key, value in input_dict.items():
    # make array for 
    y_array = len(input_dict[key])*[1]
    
    if key.endswith("/A") and range_dict[key][0] > 1:
        # mid = 10**range_dict[key][0]
        minn = range_dict[key][0]
        maxx = range_dict[key][1]
        
        # plt.scatter(mid, 1)
        plt.scatter([minn, maxx], [1,1])
        plt.scatter(input_dict[key], y_array, s = 10)
        
    elif key.endswith("/E0") or key.endswith("Ea"):
        ea2ev = 9.6e4                              
        # mid = range_dict[key][0]/ea2ev
        minn = range_dict[key][0]/ea2ev
        maxx = range_dict[key][1]/ea2ev
        
        plt.scatter([minn, maxx], [1,1])
        plt.scatter(input_dict[key]/ea2ev, y_array, s = 10)
                                     
    else: 
        plt.scatter(range_dict[key], [1,1])
        plt.scatter(input_dict[key], y_array, s = 10)
          
    plt.suptitle(key)
    
    pp.savefig(plt.gcf())
    plt.clf() 
    
    print(range_dict[key])







