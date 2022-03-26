import cantera as ct 
import numpy as np 
import pandas as pd 


class scr: 
    '''
    Single crystal facet batch reactor, modelled after the one proposed by
    C. T. Cambpbell (DOI: 10.1016/S0360-0564(08)60016-4)

    A small (~30 mL) ultrahigh vacuum chamber is prepared with a sample of the
    metal facet (~1 cm^2) to be tested. the reactant gases are introduced to the reactor 
    and the metal is brought up to temp. After 5% conversion is achieved for a
    reactant (specified in the experiment, e.g. CO methanation mechanism will 
    terminate after 5% by mol CO is consumed).


    '''
    def __init__(
        self,
        T,
        P,
        species,
        V,
        outputs,
        ):
        '''
        define the reactor 
        catalyst_area - area of the catalyst surface exposed (m^2)

        T - temperature in K
        P - pressure in Pa
        species - (dict) the species involved in the system and their
        respective mol fractions.
        V - the reactor volume in m^2
        outputs - (dict) TOFs or other benchmarks for output species 
        '''
        self.T = T
        self.P = P
        self.species = species
        self.V = V
        self.outputs = outputs


        

        















# load initial data from yang 2010
t_mult = 1.0
p_mult = 1e0
temp = [525, 550, 575, 600]
modified_temps = [
    525*t_mult, 
    550*t_mult, 
    575*t_mult, 
    600*t_mult
]
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
meoh_rates = dict(zip(temp, meoh_rates_cm))

rwgs_rates_cm = np.exp(rwgs_ln_rate)*10**15
rwgs_rates = dict(zip(temp, rwgs_rates_cm))

# get yang site density from reported rates
meoh_rate_tof = 6.3e-3
rwgs_rate_tof = 1.8
meoh_rate_cm = meoh_rates[575]
rwgs_rate_cm = rwgs_rates[575]
yang_site_density = np.mean([meoh_rate_cm/meoh_rate_tof, rwgs_rate_cm/rwgs_rate_tof])

# yang TOFs
meoh_rates_tof = {temp:rate/yang_site_density for temp,rate in meoh_rates.items()}
rwgs_rates_tof = {temp:rate/yang_site_density for temp,rate in rwgs_rates.items()}

# convert yang pressures to pascals
p_co2 = 0.5 * ct.one_atm*p_mult
p_h2 =  4.5 * ct.one_atm*p_mult
p_total = p_co2+p_h2

# get total pressure at initial CO2/H2 temp using ig law pv = nrt
initial_temp = 300 #[k]
p_total_at_temp = np.array(temp)*p_total/initial_temp

x_co2 = p_co2/p_total
x_h2 = p_h2/p_total

# create thermo phases
toy_mech = False

# specify cat area. should be somewhat arbitrary but test to be sure.
# set to 1cm^2 per campbell paper
cat_area = 1e-4


if toy_mech: 
    yaml_file = "../../External_data/park_et_al_model_reconstruction/park_mech.cti"

    gas = ct.Solution(yaml_file, "gas")
    surf = ct.Interface(yaml_file,"surface1", [gas])
    surf.coverages = {'*':1.0}

    for name in gas.species_names:
        if name.startswith("CO2"):
            co2_str = name
        if name.startswith("H2"):
            h2_str = name
        if name.startswith("CH3OH"):
            meoh_str = name
    surf_site_str = '*'
else: 
    yaml_file = "/work/westgroup/ChrisB/_01_MeOH_repos/meOH-synthesis/base/cantera/chem_annotated.cti"
    gas = ct.Solution(yaml_file, "gas")
    surf = ct.Interface(yaml_file,"surface1", [gas])

    for name in gas.species_names:
        if name.startswith("CO2("):
            co2_str = name
        if name.startswith("H2("):
            h2_str = name
        if name.startswith("CH3OH("):
            meoh_str = name
    surf_site_str = 'X(1)'


mole_fractions = {co2_str:x_co2,h2_str:x_h2}

# molecular weights for mass flow calculations
mw_co = 28.01e-3  # [kg/mol]
mw_co2 = 44.01e-3  # [kg/mol]
mw_h2 = 2.016e-3  # [kg/mol]
mw_h2o = 18.01528e-3  # [kg/mol]

# run the yang batch reactor
nspecies = len(gas.species_names)
meoh_gas_index = gas.kinetics_species_index(meoh_str)
meoh_surf_index = surf.kinetics_species_index(meoh_str)

# data from all runs 
meoh_rop_dict = {}  
meoh_moles_norm_dict = {}
temps_dict = {}
press_dict = {}
times_dict = {}
conversions_dict = {}
meoh_slope_dict = {}

reactime = 200*60
dt = 5.0

for i in range(len(temp)):
    print(temp[i])
    print("iteration", i)
    
    # data for this run
    meoh_rop = []  
    meoh_moles_norm = []
    temps = []
    press = []
    times = []
    conversions = []
    
    # get total pressure at initial CO2/H2 temp using ig law pv = nrt
    initial_temp = 300 #[k]
    p_total_at_temp = np.array(temp)*p_total/initial_temp

    x_co2 = p_co2/p_total
    x_h2 = p_h2/p_total

    # set the thermophase (maybe there is a faster way to set this up each time)
#     gas.TPX = temp[i], p_total_at_temp[i], mole_fractions
    gas.TPX = 300, p_total_at_temp[i], mole_fractions
    surf.TP = temp[i], p_total_at_temp[i]
    surf.coverages = {surf_site_str:1.0}
    r = ct.IdealGasReactor(gas, energy="off")
    r.volume = 3e-5 # from cambell at al paper, 30 mL microreactor? 
    rsurf = ct.ReactorSurface(surf, r, A=cat_area)

    # initialize reactor network
    sim = ct.ReactorNet([r])
    sim.atol = 1e-15
    sim.rtol = 1e-13

    # see if the Rop is constant
    conversion = 0
    t = 0
    while t < reactime and conversion < 0.05:
        t += dt
        sim.advance(t)
        meoh_total_rop = gas.net_production_rates[meoh_gas_index] + surf.net_production_rates[meoh_surf_index]
        meoh_total_rop = meoh_total_rop/surf.site_density
        meoh_rop.append(meoh_total_rop)
        times.append(sim.time/60) # time in minutes
        temps.append(gas.T)
        press.append(gas.P)

        moles_co2_initial = p_co2*r.volume/(ct.gas_constant*initial_temp) # Pco2*v/RT
        moles_meoh = (gas[meoh_str].X*gas.P*r.volume)/(ct.gas_constant*gas.T) # mole frac*total moles
        moles_co2 = (gas[co2_str].X*gas.P*r.volume)/(ct.gas_constant*gas.T) # mole frac*total moles
        
        
        # calculate conversion:
        # (moles meoh possible (starting moles co2)-moles meoh current)/moles co2 initial
        conversion = (moles_co2_initial - moles_co2)/moles_co2_initial
        conversions.append(conversion)

        # calculate the moles methanol normalized
        # to the number of surface sites
        meoh_moles_norm.append(float(moles_meoh/(surf.site_density*cat_area))) 

    meoh_rop_dict[temp[i]] = meoh_rop
    meoh_moles_norm_dict[temp[i]] = meoh_moles_norm
    temps_dict[temp[i]] = temps
    press_dict[temp[i]] = press
    times_dict[temp[i]] = times
    conversions_dict[temp[i]] = conversions

    print(len(times), len(meoh_moles_norm))
    meoh_slope, intercept, r_value, p_value, std_err = stats.linregress(times,meoh_moles_norm)
    meoh_slope_dict[temp[i]] = meoh_slope
    print(
        "slope: ", meoh_slope,
        "\nintercept: ", intercept, 
        "\nr^2 value: ", r_value**2, 
        "\np value: ", p_value, 
        "\nstd error: ", std_err 
    )
    plt.plot(times,conversions)
    plt.show()

    
# build_and_run_yang()