import yaml
import numpy as np
import pandas as pd
import cantera as ct
import itertools
import scipy.io
import math
import sys
import os
from pathlib import Path

repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# manually give addresses for data.
# graaf_data_dir = '/work/westgroup/ChrisB/_01_MeOH_repos/meOH-analysis/cantera_simulations/Graaf_data/'
# yang_data_dir = '/work/westgroup/ChrisB/_01_MeOH_repos/meOH-analysis/cantera_simulations/yang_2010_data/'
# grabow_conditions_dir = '/work/westgroup/ChrisB/_01_MeOH_repos/meOH-analysis/cantera_simulations/Grabow_data/original_runs/'

# # make addressed arguments
# graaf_data_dir = sys.argv[1]
# # yang_data_dir = sys.argv[2]
# grabow_conditions_dir = sys.argv[2]

graaf_data_dir = os.path.join(repo_dir, "experimental_data/graaf")
grabow_conditions_dir = os.path.join(repo_dir, "experimental_data/grabow")

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

###############################################################################
# Grabow Spinning Basket reactor data
# pull in data from Grabow's study of the effect of CO2/CO/H2/H2O concentration
# on the output of 
###############################################################################

# make conditions dict for grabow runs to start: 
conditions_dict = {}

# get vol
volume = ((35e-3)**2)*np.pi*(70e-3)/2

# catalyst area
site_density = 5*61.67*1e-6*1e3 # [moles/kg]

total_sites = site_density*4.24e-3 #moles sites (4.24 g cat)

rmg_site_density_cu = 2.943e-9*1e4 #mol/m^2 see chemkin surface file

cat_area = (total_sites)/(rmg_site_density_cu) #mol/mol/m^2()

# Volume flow rate
Vin_cm3_min = 470.4
Vin_m3_sec = Vin_cm3_min/(60*10**6)
volume_flow = Vin_m3_sec

temp = 528 # K
pressure = 75*101325 # atm
expt_name = 'grabow2011'
expt_type = 'sbr'

# mole fractions 
CO2_ratio = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] #CO2/(CO+CO2)
H2_moles = [0.5, 0.75, 0.8, 0.95]
H2O_moles = [0.05]

#get total number of runs
size = len(H2_moles)*len(CO2_ratio)

H2_mole_list = []
CO2_mole_list = []
CO_mole_list = []
H2O_mole_list = []

# for h2 in H2_moles: 
#     for co2_r in CO2_ratio: 
#         if h2 == 0.75:
#             h2o = 0.05
#         else: 
#             h2o = 0.
        
#         co2 = (1 - h2o - h2)*co2_r
#         co2 = round(co2,3)
#         co = (1 - h2o - h2)*(1-co2_r)
#         co = round(co, 3)
#         H2_mole_list.append(h2)
#         CO2_mole_list.append(co2)
#         CO_mole_list.append(co)
#         H2O_mole_list.append(h2o)
            

# mole_dict = {
#     'H2':H2_mole_list,
#     'CO2':CO2_mole_list,
#     'CO':CO_mole_list,
#     'H2O':H2O_mole_list,
# }

# # make mole dict into a list of dictionaries for easier implementation
# mole_dict = list(zip_dict(**mole_dict))

# # check that mole fractions add to 1
# for i in range(len(H2_mole_list)):
#     summy = mole_dict[i]['H2'] + mole_dict[i]['CO2']+mole_dict[i]['CO']+mole_dict[i]['H2O']
#     if summy !=1:
#         print(f"mole fractions do not add up to one! add up to {summy}") 

# make an output dict 
# for grabow, currently only care about the 
# meoh TOF and the H2O TOF

Grabow_rates = pd.read_csv(os.path.join(grabow_conditions_dir, "Grabow_rates.csv"))
grabow_input_labels = {
    "y(CO)":"CO",
    "y(CO2)":"CO2",
    "y(H2)":"H2",
    "y(H2O)":"H2O",
}

mole_dict = {}
for df_key,inp in grabow_input_labels.items():
    mole_dict[inp] = list(Grabow_rates[df_key])

# make mole dict into a list of dictionaries for easier implementation
mole_dict = list(zip_dict(**mole_dict))

# check that mole fractions add to 1
for i in range(len(H2_mole_list)):
    summy = mole_dict[i]['H2'] + mole_dict[i]['CO2']+mole_dict[i]['CO']+mole_dict[i]['H2O']
    if summy !=1:
        print(f"mole fractions do not add up to one! add up to {summy}")
        print(mole_dict[i]['H2'], mole_dict[i]['CO2'], mole_dict[i]['CO'], mole_dict[i]['H2O'])  
        
# get outputs 
grabow_output_labels = {
    "Methanol Production":"CH3OH",
    "H2O Production":"H2O",
}

output_dict = {}
for df_key,inp in grabow_output_labels.items():
    output_dict[inp] = list(Grabow_rates[df_key])
    
output_dict = list(zip_dict(**output_dict))

# make yaml
grabow_yammy = {}
grabow_yammy['expt_name'] = [expt_name]*size
grabow_yammy['volume'] = [volume]*size
grabow_yammy['catalyst_area']= [cat_area]*size
grabow_yammy['volume_flowrate'] = [volume_flow]*size
grabow_yammy['temperature'] = [temp]*size
grabow_yammy['pressure'] = [pressure]*size
grabow_yammy['experiment_type'] = [expt_type]*size
grabow_yammy['species'] = mole_dict
grabow_yammy['output'] = output_dict

print("temperature dtype ",  type(temp))

# make determination of error under 40% 
# load in Graaf experimental data
# exclude feed 8 because it is a monolith reactor
path_str = os.path.join(graaf_data_dir,"combined_experimental_runs.xlsx")
df_graaf = pd.read_excel(path_str, engine='openpyxl')


# load in initial conditions from experiment from matlab
translation = {39: None, 91: None, 93: None}  # remove ', [, and ]
count = 0
first_row = True

for path in Path(grabow_conditions_dir).rglob('*.mat'):
    # get feed number and run number from file name
    path_str = str(path)
    # need a way to match the grabow experiments with the graaf ones. their numbering scheme is nonsense.
    file_name = path_str.replace(".mat", "").split("/")[-1]
    feed = int(file_name.split("_")[2])
#     run = int(file_name.split("_")[4].split(".")[0])

    mat = scipy.io.loadmat(path_str)
    conditions = mat['ccondition']
    reactions = mat['reaction']

    # get relevant TOFs, catalyst weights, # sites
    meoh_TOF = float(conditions["MeOH_TOF"])
    h2o_TOF = float(conditions["H2O_TOF"])
    cat_weight = float(conditions["CatalystWeight"])
    num_sites = float(conditions["Sites"])
    # molar flow rate divided by moles sites
    flow_site_basis = float(conditions["Fin"])
    flow_cm_3_min = float(conditions["Flow"])
    P_atm = float(conditions["P"])
    T_grabow_K = float(conditions["T"])
    inlet_CO = float(conditions["molfrac"][0][0][0][0][1][0][0])
    inlet_CO2 = float(conditions["molfrac"][0][0][0][0][1][0][1])
    inlet_h2 = float(conditions["molfrac"][0][0][0][0][1][0][2])

    #read in species names to use as columns in dataframe
    species_names = []

    for spec in range(len(mat['species'][0, :])):
        species_string = str(mat['species'][0, spec][0])
        species_string = species_string.translate(translation)
        species_names.append(species_string)

    # read in results (surface coverages, partial pressures)
    results = mat['Y']

    if first_row:
        # create data frame with species names as column headers
        df_grabow = pd.DataFrame(data=results, columns=species_names)
        df_grabow = df_grabow.tail(1)
        df_grabow["feed"] = feed
#         df_grabow["run"] = run
        df_grabow["MeOH TOF (1/s)"] = meoh_TOF
        df_grabow["H2O TOF (1/s)"] = h2o_TOF
        df_grabow["Catalyst Weight (g)"] = cat_weight
        df_grabow["Number of Sites (mol)"] = num_sites
        df_grabow["Fin site basis (mol/site/s)"] = flow_site_basis
        df_grabow["flow (cm^3/min)"] = flow_cm_3_min
        df_grabow["P (atm)"] = P_atm
        df_grabow["T grab (K)"] = T_grabow_K
        df_grabow["inlet mol frac CO"] = inlet_CO
        df_grabow["inlet mol frac CO2"] = inlet_CO2
        df_grabow["inlet mol frac h2"] = inlet_h2

        # get total pressure by adding each species partial pressure
        total_pressure = 0
        for column in df_grabow:
            if "g" in column[-1].strip():
                total_pressure += float(df_grabow[column])
        df_grabow["total pressure (bar)"] = total_pressure

        first_row = False
    else:
        # add tail from new dataframe
        new_df = pd.DataFrame(data=results, columns=species_names)
        new_df = new_df.tail(1)
        new_df["feed"] = feed
#         new_df["run"] = run
        new_df["MeOH TOF (1/s)"] = meoh_TOF
        new_df["H2O TOF (1/s)"] = h2o_TOF
        new_df["Catalyst Weight (g)"] = cat_weight
        new_df["Number of Sites (mol)"] = num_sites
        new_df["Fin site basis (mol/site/s)"] = flow_site_basis
        new_df["flow (cm^3/min)"] = flow_cm_3_min
        new_df["P (atm)"] = P_atm
        new_df["T grab (K)"] = T_grabow_K
        new_df["inlet mol frac CO"] = inlet_CO
        new_df["inlet mol frac CO2"] = inlet_CO2
        new_df["inlet mol frac h2"] = inlet_h2

        # get total pressure by adding each species partial pressure
        total_pressure = 0
        for column in new_df:
            if "g" in column[-1].strip():
                total_pressure += float(new_df[column])

        new_df["total pressure (bar)"] = total_pressure

        df_grabow = df_grabow.append(new_df, ignore_index=True)


# convert all partial pressures to mole fractions
# any value with a "g" (gas) after the species name
for column in df_grabow:
    if "g" in column[-1].strip():
        df_grabow[column] = df_grabow[column]/df_grabow["total pressure (bar)"]

# calculate turn over frequencies
df_grabow["Grabow MeOH TOF (1/s)"] = df_grabow["CH3OHg"] * \
    df_grabow["Fin site basis (mol/site/s)"]
df_grabow["H2O TOF (1/s)"] = df_grabow["H2Og"] * \
    df_grabow["Fin site basis (mol/site/s)"]

count = 0
found = False
for indexgrabow, row_grabow in df_grabow.iterrows():
    for index_graaf, row_graaf in df_graaf.iterrows():

        # get the values in the Graaf experiment that match the grabow values to 1e-3
        if (math.isclose(row_graaf["10^6 * V (M^3/s)"], row_grabow["flow (cm^3/min)"]/60, rel_tol=1e-3) and
            math.isclose(row_graaf["T(K)"], row_grabow["T grab (K)"], rel_tol=1e-3) and
            math.isclose(row_graaf["p (bar)"], row_grabow["P (atm)"]*1.01325, rel_tol=1e-3) and
            math.isclose(row_graaf["feed Yco"], row_grabow["inlet mol frac CO"], rel_tol=1e-3) and
            math.isclose(row_graaf["feed Yco2"], row_grabow["inlet mol frac CO2"], rel_tol=1e-3) and
            math.isclose(row_graaf["feed Yh2"], row_grabow["inlet mol frac h2"], rel_tol=1e-3) and
                math.isclose(row_graaf["MeOH TOF (mol/site/s)"], row_grabow["MeOH TOF (1/s)"], rel_tol=1e-3)):

            new_series_data = row_grabow.append(row_graaf)
            found = True
    if found:
        if count == 0:
            df_combined = pd.DataFrame(
                new_series_data, index=new_series_data.index).transpose()
            count += 1

        if count > 0:
            df_combined = df_combined.append(
                pd.DataFrame(new_series_data).transpose())
            count += 1

        found = False

# translation table
trans_table = {}
error_under_forty = list(df_combined['Grabow_exp_index'])


###############################################################################
# Graaf Data experimental 
###############################################################################

file_name_feed1 = os.path.join(graaf_data_dir,"Feed_1.xlsx")
file_name_feed2 = os.path.join(graaf_data_dir,"Feed_2.xlsx")
file_name_feed3 = os.path.join(graaf_data_dir,"Feed_3.xlsx")
file_name_feed4 = os.path.join(graaf_data_dir,"Feed_4.xlsx")
file_name_feed5 = os.path.join(graaf_data_dir,"Feed_5.xlsx")
file_name_feed6a = os.path.join(graaf_data_dir,"Feed_6a.xlsx")
file_name_feed6b = os.path.join(graaf_data_dir,"Feed_6b.xlsx")
file_name_feed7a = os.path.join(graaf_data_dir,"Feed_7a.xlsx")
file_name_feed7b = os.path.join(graaf_data_dir,"Feed_7b.xlsx")
file_name_all = os.path.join(graaf_data_dir,"combined_experimental_runs.xlsx")

df_1 = pd.read_excel(file_name_feed1, engine='openpyxl')
df_2 = pd.read_excel(file_name_feed2, engine='openpyxl')
df_3 = pd.read_excel(file_name_feed3, engine='openpyxl')
df_4 = pd.read_excel(file_name_feed4, engine='openpyxl')
df_5 = pd.read_excel(file_name_feed5, engine='openpyxl')
df_6a = pd.read_excel(file_name_feed6a, engine='openpyxl')
df_6b = pd.read_excel(file_name_feed6b, engine='openpyxl')
df_7a = pd.read_excel(file_name_feed7a, engine='openpyxl')
df_7b = pd.read_excel(file_name_feed7b, engine='openpyxl')
df = pd.read_excel(file_name_all, engine='openpyxl')

# Needed: [T, P, V, YH2, YCO2, wcat] -- Create a list of lists
# Should be columns 2 (T), 1 (P), 3 (V), 6 (YH2), 5 (YCO2) 6 (cat weight)
# Each list is the conditions of one experimental Graaf run

# List of dataframes with feed conditions
df_list = [df_1, df_2, df_3, df_4, df_5, df_6a, df_6b, df_7a, df_7b]

# Loop through dataframes and create a list of conditions based on Graaf runs
# Loop through each row in the dataframes and add that row's conditions to the list of lists

H2_mole_list = []
CO2_mole_list = []
CO_mole_list = []

H2_mole_list_out = []
CO2_mole_list_out = []
CO_mole_list_out = []
CH3OH_mole_list_out = []
H2O_mole_list_out = []

volume_flows = []
cat_areas = []
pressures = []
temps = []
meoh_tofs = []
h2o_tofs = []
run_nums = []
use_for_opt = []

# for i in range(len(df_list)):
#     df = df_list[i]
for row in range(len(df)):
    if not np.isnan(df.iloc[row, df.columns.get_loc('T(K)')]):
        
        # moles
        H2_mole_list.append(float(df.iloc[row,df.columns.get_loc('feed Yh2')]))
        CO2_mole_list.append(float(df.iloc[row,df.columns.get_loc('feed Yco2')]))
        CO_mole_list.append(float(df.iloc[row,df.columns.get_loc('feed Yco')]))

        # moles out
        H2_mole_list_out.append(float(df.iloc[row,df.columns.get_loc('Yh2')]))
        CO2_mole_list_out.append(float(df.iloc[row,df.columns.get_loc('Yco2')]))
        CO_mole_list_out.append(float(df.iloc[row,df.columns.get_loc('Yco')]))
        CH3OH_mole_list_out.append(float(df.iloc[row,df.columns.get_loc('Ych3oh')]))
        H2O_mole_list_out.append(float(df.iloc[row,df.columns.get_loc('Yh2o')]))
        
        # volume flow
        volume_flow = float(df.iloc[row,df.columns.get_loc('10^6 * V (M^3/s)')])*1e-6  # m^3
        volume_flows.append(volume_flow)   
        
        # catalyst weight
        cat_weight = float(df.iloc[row,df.columns.get_loc('wcat (g)')])*1e-3 # [kg]
        cat_area = (cat_weight * site_density)/rmg_site_density_cu  # [m^3]
        cat_areas.append(cat_area)
        
        # Pressure
        pressure = float(df.iloc[row,df.columns.get_loc('p (bar)')])*1e5 #[Pa]
        pressures.append(pressure)
        
        # temperatures
        temps.append(float(df.iloc[row,df.columns.get_loc('T(K)')]))  
        
        # meoh TOF
        meoh_tofs.append(float(df.iloc[row,df.columns.get_loc('MeOH TOF (mol/site/s)')]))
        
        # h2o TOF
        h2o_tofs.append(float(df.iloc[row,df.columns.get_loc('H2O TOF (mol/site/s)')]))


        # feed # (used to determine if it was used by grabow)
        run_num = float(df.iloc[row,df.columns.get_loc('Grabow_exp_index')])
        run_nums.append(run_num)
        
        if run_num in error_under_forty:
            use_for_opt.append(True)
        else: 
            use_for_opt.append(False)
                                        
# check that mole fractions add to 1
for i in range(len(H2_mole_list)):
    summy = H2_mole_list[i] + CO2_mole_list[i]+CO_mole_list[i]
    if summy !=1:
        print(f"mole fractions do not add up to one! add up to {summy}")              

size = len(temps)

# construct mole dict
mole_dict = {
    'H2':H2_mole_list,
    'CO2':CO2_mole_list,
    'CO':CO_mole_list,
    'H2O':[0]*size
}  

# make mole dict into a list of dictionaries for easier implementation
mole_dict = list(zip_dict(**mole_dict))

# construct output mole_dict
mole_out_dict = {
    'H2':H2_mole_list_out,
    'CO2':CO2_mole_list_out,
    'CO':CO_mole_list_out,
    'CH3OH':CH3OH_mole_list_out,
    'H2O':H2O_mole_list_out,
}

mole_out_dict = list(zip_dict(**mole_out_dict))

output_dict = {
    'CH3OH':meoh_tofs,
    'H2O':h2o_tofs,
}

output_dict = list(zip_dict(**output_dict))

expt_name = "graaf_1988"
# make yaml
graaf_yammy = {}
graaf_yammy['expt_name'] = [expt_name]*size
graaf_yammy['volume'] = [volume]*size
graaf_yammy['catalyst_area']= cat_areas
graaf_yammy['volume_flowrate'] = volume_flows
graaf_yammy['temperature'] = temps
graaf_yammy['pressure'] = pressures
graaf_yammy['experiment_type'] = [expt_type]*size
graaf_yammy['species'] = mole_dict  
graaf_yammy['species_out'] = mole_out_dict
graaf_yammy['output'] = output_dict  
graaf_yammy['use_for_opt'] = use_for_opt


###############################################################################
# yang batch reactor experiments
###############################################################################

# Load in values from plot in Yang 2010
temps = [525, 550, 575, 600]

size = len(temps)

meoh_ln_rate = [
    -6.691144708,
    -5.978401728,
    -4.48812095,
    -3.894168467,
]

rwgs_ln_rate = [
    -0.578342066,
    0.572607525,
    1.171517945,
    2.072487534,
]

# convert to molecules/cm^2/sec
meoh_rates_cm = np.exp(meoh_ln_rate)*10**15
meoh_rates = dict(zip(temps, meoh_rates_cm))

rwgs_rates_cm = np.exp(rwgs_ln_rate)*10**15
rwgs_rates = dict(zip(temps, rwgs_rates_cm))

meoh_rate_tof = 6.3e-3
rwgs_rate_tof = 1.8
meoh_rate_cm = meoh_rates[575]
rwgs_rate_cm = rwgs_rates[575]
site_density = np.mean([meoh_rate_cm/meoh_rate_tof, rwgs_rate_cm/rwgs_rate_tof])
"site density molecules/cm: {:2e}".format(site_density)

meoh_rates_tof = {temp:rate/site_density for temp,rate in meoh_rates.items()}
rwgs_rates_tof = {temp:rate/site_density for temp,rate in rwgs_rates.items()}

# recast values to floats, because yaml cannot handle objects for some reason
for temp in meoh_rates_tof.keys():
    meoh_rates_tof[temp] = float(meoh_rates_tof[temp])
    rwgs_rates_tof[temp] = float(rwgs_rates_tof[temp])

# convert to pascals
p_co2 = 0.5 * ct.one_atm
p_h2 =  4.5 * ct.one_atm
p_total = p_co2+p_h2

# get total pressure at temp using ig law pv = nrt
initial_temp = 300 #[k]
p_total_at_temp = np.array(temps)*p_total/initial_temp

# get mole fractions
x_co2 = p_co2/p_total
x_h2 = p_h2/p_total


H2_mole_list = [x_h2]*size
CO2_mole_list = [x_co2]*size

mole_dict = {
    'H2':H2_mole_list,
    'CO2':CO2_mole_list,
}  

# make mole dict into a list of dictionaries for 
# easier implementation
mole_dict = list(zip_dict(**mole_dict))

output_dict = {
    'CH3OH':list(meoh_rates_tof.values()),
    'CO':list(rwgs_rates_tof.values()),
}  

# make mole dict into a list of dictionaries for 
# easier implementation
output_dict = list(zip_dict(**output_dict))
expt_name = "yang_2010"

# make yaml
yang_yammy = {}
yang_yammy['expt_name'] = [expt_name]*size
yang_yammy['volume'] = [0.1]*size
yang_yammy['catalyst_area']= [1e2]*size
yang_yammy['temperature'] = temps
yang_yammy['pressure'] = p_total_at_temp.tolist()
yang_yammy['experiment_type'] = ["batch"]*size
yang_yammy['species'] = mole_dict
yang_yammy['output'] = output_dict

yaml_output = {}
yaml_output['yang_2010'] = yang_yammy
yaml_output['graaf_1988'] = graaf_yammy
yaml_output['grabow_2011'] = grabow_yammy

output_file = "all_experiments.yaml"
with open(output_file , 'w') as f:
    doc = yaml.safe_dump(yaml_output, f)
