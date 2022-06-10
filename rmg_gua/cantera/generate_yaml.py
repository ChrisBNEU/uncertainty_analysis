import yaml
import numpy as np
import pandas as pd
import cantera as ct
import itertools

# manually give addresses for data.
graaf_data_dir = '/work/westgroup/ChrisB/_01_MeOH_repos/meOH-analysis/cantera_simulations/Graaf_data/'
yang_data_dir = '/work/westgroup/ChrisB/_01_MeOH_repos/meOH-analysis/cantera_simulations/yang_2010_data/'
grabow_conditions_dir = '/work/westgroup/ChrisB/_01_MeOH_repos/meOH-analysis/cantera_simulations/Grabow_data/original_runs/'

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

Grabow_rates = pd.read_csv("/work/westgroup/ChrisB/_01_MeOH_repos/meOH-analysis/cantera_simulations/Grabow_data/Paper_plot_data/Grabow_rates.csv")
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
###############################################################################
# Graaf Data experimental 
###############################################################################

file_name_feed1 = graaf_data_dir + "Feed_1.xlsx"
file_name_feed2 = graaf_data_dir + "Feed_2.xlsx"
file_name_feed3 = graaf_data_dir + "Feed_3.xlsx"
file_name_feed4 = graaf_data_dir + "Feed_4.xlsx"
file_name_feed5 = graaf_data_dir + "Feed_5.xlsx"
file_name_feed6a = graaf_data_dir + "Feed_6a.xlsx"
file_name_feed6b = graaf_data_dir + "Feed_6b.xlsx"
file_name_feed7a = graaf_data_dir + "Feed_7a.xlsx"
file_name_feed7b = graaf_data_dir + "Feed_7b.xlsx"

df_1 = pd.read_excel(file_name_feed1, engine='openpyxl')
df_2 = pd.read_excel(file_name_feed2, engine='openpyxl')
df_3 = pd.read_excel(file_name_feed3, engine='openpyxl')
df_4 = pd.read_excel(file_name_feed4, engine='openpyxl')
df_5 = pd.read_excel(file_name_feed5, engine='openpyxl')
df_6a = pd.read_excel(file_name_feed6a, engine='openpyxl')
df_6b = pd.read_excel(file_name_feed6b, engine='openpyxl')
df_7a = pd.read_excel(file_name_feed7a, engine='openpyxl')
df_7b = pd.read_excel(file_name_feed7b, engine='openpyxl')

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
volume_flows = []
cat_areas = []
pressures = []
temps = []
meoh_tofs = []
h2o_tofs = []

for i in range(len(df_list)):
    df = df_list[i]
    for row in range(len(df)):
        if not np.isnan(df.iloc[row, df.columns.get_loc('T(K)')]):
            
            # moles
            H2_mole_list.append(float(df.iloc[row,df.columns.get_loc('feed Yh2')]))
            CO2_mole_list.append(float(df.iloc[row,df.columns.get_loc('feed Yco2')]))
            CO_mole_list.append(float(df.iloc[row,df.columns.get_loc('feed Yco')]))
            
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
graaf_yammy['output'] = output_dict  


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
