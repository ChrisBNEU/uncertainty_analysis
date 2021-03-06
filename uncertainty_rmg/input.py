# Running on current commits 13Jul2021 for 
# RMG-Py: forbidden_input
# RMG-Database: meoh_3
# running with david's forbidden input branch for Py, 
# Database:
# bjarne's Pt111 thermo
# the correct H2vdw values for thermo per david's comments on PR #516 
# removed grabow rates from training data 

# Data sources
database(
    thermoLibraries=[
        # 'surfaceThermoCu111', 
        'surfaceThermoPt111', 
        'primaryThermoLibrary', 
        'thermo_DFT_CCSDTF12_BAC',
        'DFT_QCI_thermo'
        ],
    reactionLibraries = [
        'BurkeH2O2inArHe',
        'BurkeH2O2inN2',
        # comment out libraries until we can perturb them too. 
        # 'Surface/Methane/Deutschmann_Ni',
        'Surface/CPOX_Pt/Deutschmann2006_adjusted',
        ],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies =[('surface',True),'default'],
    kineticsEstimator = 'rate rules',
)

catalystProperties( # default values for Cu(111) calculated by Katrin Blondal and Bjarne Kreitz at Brown University
    bindingEnergies = {
                       'C':(-4.96033553, 'eV/molecule'),
                       'O':(-4.20763879, 'eV/molecule'),
                       'N':(-3.58446699, 'eV/molecule'),
                       'H':(-2.58383235, 'eV/molecule'),
                       },
    surfaceSiteDensity=(2.943e-9, 'mol/cm^2'),  # from Katrin
    coverageDependence=True,

)


# List of species
species(
    label='X',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N"),
)

species(
    label='H2',
    reactive=True,
    structure=SMILES("[H][H]"),
)

species(
    label='CO',
    reactive=True,
    structure=SMILES("[C-]#[O+]"),
)

species(
    label='CO2',
    reactive=True,
    structure=SMILES("O=C=O"),
)

species(
    label='H2O',
    reactive=True,
    structure=SMILES("O"),
)

species(
    label='CH3OH',
    reactive=True,
    structure=SMILES("CO"),
)

#----------
# Reaction systems
surfaceReactor(
    temperature=[(400,'K'),(700, 'K')],
    initialPressure=(15.0, 'bar'),
    nSims = 4,
    initialGasMoleFractions={
        "CO": 0.0,
        "CO2": 0.776,
        "H2": 0.7669,
        "N2": 0.1555,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO2":0.99,},
    terminationTime=(10., 's'),
    terminationRateRatio=0.01
)

surfaceReactor(
    temperature=[(400,'K'),(700, 'K')],
    initialPressure=(76.0, 'bar'),
    nSims = 4,
    initialGasMoleFractions={
        "CO": 0.0,
        "CO2": 0.776,
        "H2": 0.7669,
        "N2": 0.1555,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO2":0.99,},
    terminationTime=(10., 's'),
    terminationRateRatio=0.01
)

surfaceReactor(
    temperature=[(400,'K'),(700, 'K')],
    initialPressure=(15.0, 'bar'),
    nSims = 4,
    initialGasMoleFractions={
        "CO": 0.2019,
        "CO2": 0.0,
        "H2": 0.6424,
        "N2": 0.1557,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO":0.99,},
    terminationTime=(10., 's'),
    terminationRateRatio=0.01
)

surfaceReactor(
    temperature=[(400,'K'),(700, 'K')],
    initialPressure=(76.0, 'bar'),
    nSims = 4,
    initialGasMoleFractions={
        "CO": 0.2019,
        "CO2": 0.0,
        "H2": 0.6424,
        "N2": 0.1557,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO":0.99,},
    terminationTime=(10., 's'),
    terminationRateRatio=0.01
)

# adding in reactor at 15 and 75 bar with all starting species present because why not. 
surfaceReactor(
    temperature=[(400,'K'),(700, 'K')],
    initialPressure=(15.0, 'bar'),
    nSims = 4,
    initialGasMoleFractions={
        "CO": 0.2,
        "CO2": 0.2,
        "H2": 0.2,
        "H2O": 0.2,
        "N2": 0.1,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO":0.99,},
    terminationTime=(10., 's'),
    terminationRateRatio=0.01
)

surfaceReactor(
    temperature=[(400,'K'),(700, 'K')],
    initialPressure=(76.0, 'bar'),
    nSims = 4,
    initialGasMoleFractions={
        "CO": 0.2,
        "CO2": 0.2,
        "H2": 0.2,
        "H2O": 0.2,
        "N2": 0.1,
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO":0.99,},
    terminationTime=(10., 's'),
    terminationRateRatio=0.01
)

simulator(
    atol=1e-18,
    rtol=1e-12,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.01,
# inturrupt tolerance was 0.1 wout pruning, 1e8 w pruning on
    toleranceInterruptSimulation=0.01,
    maximumEdgeSpecies=500000,
# PRUNING: uncomment to prune
#    minCoreSizeForPrune=50,
# prune before simulation based on thermo
#    toleranceThermoKeepSpeciesInEdge=0.5,
# prune rxns from edge that dont move into core
#    minSpeciesExistIterationsForPrune=2,
# FILTERING: set so threshold is slightly larger than max rate constants
#    filterReactions=True,
#    filterThreshold=5e8, # default value
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=False,
    generatePlots=False,
    saveEdgeSpecies=False,
    saveSimulationProfiles=False,
    generateSeedEachIteration=False
)

generatedSpeciesConstraints(
    allowed=['input species','reaction libraries'],
#    maximumRadicalElectrons=2,
    maximumCarbonAtoms=1,
)

forbidden(
    label='CO2_bidentate',
    structure=adjacencyList(
        """
        1 O u0 p2 c0 {2,D}
        2 C u0 p0 c0 {1,D} {3,S} {4,S}
        3 X u0 p0 c0 {2,S}
        4 O u0 p2 c0 {2,S} {5,S}
        5 X u0 p0 c0 {4,S}
        """
    )
)

# forbidden(
#     label='H2X',
#     structure=adjacencyList(
#         """
#         1 H u0 p0 c0 {2,S}
#         2 H u0 p0 c0 {1,S}
#         3 X u0 p0 c0
#         """
#     )
# )

