
import PEUQSE
import numpy as np 
import pandas as pd 
import yaml
import math
import rms_simulation



def initial_val_unc(kinetics, num): 
    """ 
    return the initial value, an estimate of parameter uncertainty, and upper
     as a tuple
    A-factors: if not sticking coeffient, +/- 1, return LogA. else, +/- 0.5
    Ea: +/- 0.3 eV. 

    Limits: 
    A (upper):  - 1 for sticking coefficient, 1e25 for others
    A (lower):  - 0 for sticking coefficient, 0 for others
    E (upper):  - 400 kJ/mol
    E (lower):  - 0 kj/mol
    """ 
    # logA if it is not sticking coefficient
    if kinetics["A"] > 1: 
        A = math.log10(kinetics["A"])
        A_unc = 1
        A_lb = 0
        A_ub = 1e25
        label_A = f"A_log_{num}"
    elif kinetics["A"] <= 1: 
        A = kinetics["A"]
        A_unc = 1
        A_lb = 0
        A_ub = 1
        label_A = f"A_stick_{num}"
    
    E = kinetics["Ea"]
    E_unc = 30000
    E_lb = 0
    E_ub = 400000
    label_E = f"E_{num}"
    
    return (A, E), (A_unc, E_unc), (A_lb, E_lb), (A_ub, E_ub), (label_A, label_E)


# get the input parameter list. we need: 
# for A, Ea, etc: 
# 1. a list of the initial values
# 2. a list of the uncertainties

# for expt data
# 3. x_data: the pressure, temp, etc for each experiment
# 4. uncertainty data for the x_data
# 4. y_data: the output mole 
# fractions for CO, CO2, etc. 
# 6. Uncertainty data for the y_data
rms_file = "../baseline/rms/chem53.rms"
with open(rms_file, "r") as f: 
    rms_mech = yaml.load(f, Loader=yaml.FullLoader)

initial_values = []
uncertainties = []
lower_bounds = []
upper_bounds = []
labels = []

for num, rxn in enumerate(rms_mech["Reactions"]):
    kin = rxn["kinetics"]
    new_kin = initial_val_unc(kin, num)
    initial_values.extend(new_kin[0])
    uncertainties.extend(new_kin[1])
    lower_bounds.extend(new_kin[2])
    upper_bounds.extend(new_kin[3])
    labels.extend(new_kin[4])


# make a new yaml for these lists 
yaml_dict = {
    "initial_values": initial_values,
    "uncertainties": uncertainties,
    "lower_bounds": lower_bounds,
    "upper_bounds": upper_bounds,
    "labels": labels
}
with open("rms_initial.yaml", "w") as f: 
    yaml.dump(yaml_dict, f)

print(len(yaml_dict["initial_values"]))



# now get experimental data
expt_yaml = "../gua_cantera/all_experiments_reorg_sbr.yaml"

# get all experiments that have "use_for_opt" = True and Temp < 518L

with open(expt_yaml, "r") as f:
    expt_list = yaml.load(f, Loader=yaml.FullLoader)

new_expt_yaml = []
for expt in expt_list:
    if "use_for_opt" in expt.keys():
        if expt["use_for_opt"] == True and expt["temperature"] < 518:
            new_expt_yaml.append(expt)

# save the new yaml
new_expt_file = "./expt_data.yaml"
with open(new_expt_file, "w") as f:
    yaml.dump(new_expt_yaml, f)


# load the experimental data yaml so we get a dict of lists for each parameter
expt_data = rms_simulation.repackage_yaml(new_expt_file)
print(len(expt_data['catalyst_area']))

# save the yaml file 
with open("expt_data.yaml", "w") as f:
    yaml.dump(expt_data, f)

expt_data.keys()






