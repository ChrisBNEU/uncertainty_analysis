# minimal example RMG_GUA

# Data sources
database(
    thermoLibraries=[
        'surfaceThermoPt111', 
        'primaryThermoLibrary', 
        'thermo_DFT_CCSDTF12_BAC',
        'DFT_QCI_thermo'
        ],
    reactionLibraries = [
        'BurkeH2O2inArHe',
        'BurkeH2O2inN2',
        'Surface/CPOX_Pt/Deutschmann2006_adjusted'
        ],
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies =[('surface', True),('default', False)],
    kineticsEstimator = 'rate rules',
)

catalystProperties(
    metal = "Cu111",
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
    label='CH2O',
    reactive=True,
    structure=SMILES("C=O"),
)

species(
    label='HCOOH',
    reactive=True,
    structure=SMILES("O=CO"),
)

species(
    label='CH3OH',
    reactive=True,
    structure=SMILES("CO"),
)

species(
    label='HCOOCH3',
    reactive=True,
    structure=SMILES("O=COC"),
)

species(
   label='H*',
   reactive=True,
   structure=adjacencyList(
       """
1 H u0 p0 c0 {2,S}
2 X u0 p0 c0 {1,S}
"""),
)

species(
   label='O*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
"""),
)

species(
   label='OH*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,S}
"""),
)

species(
   label='H2O*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 X u0 p0 c0
"""),
)

species(
   label='CO*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,D}
3 X u0 p0 c0 {2,D}
"""),
)

species(
   label='CO2*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
4 X u0 p0 c0
"""),
)

species(
   label='HCO*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 X u0 p0 c0 {2,S}
"""),
)

species(
   label='HCOO*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {3,S} {5,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {4,S}
4 H u0 p0 c0 {3,S}
5 X u0 p0 c0 {1,S}
"""),
)

species(
   label='COOH*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {3,S} {4,S}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,S} {2,D} {5,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0 {3,S}
"""),
)

species(
   label='HCOOH*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {2,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {4,D}
3 H u0 p0 c0 {2,S}
4 O u0 p2 c0 {2,D}
5 X u0 p0 c0
6 H u0 p0 c0 {1,S}
"""),
)

species(
   label='CH2O*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 X u0 p0 c0
"""),
)

species(
   label='CH3O*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {2,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 X u0 p0 c0 {1,S}
"""),
)

species(
   label='CH3O2*',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {3,S} {6,S}
2 O u0 p2 c0 {3,S} {7,S}
3 C u0 p0 c0 {1,S} {2,S} {4,S} {5,S}
4 H u0 p0 c0 {3,S}
5 H u0 p0 c0 {3,S}
6 H u0 p0 c0 {1,S}
7 X u0 p0 c0 {2,S}
"""),
)

species(
   label='CH3OH*',
   reactive=True,
   structure=adjacencyList(
       """
 1 O u0 p2 c0 {2,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
7 X u0 p0 c0
"""),
)

species(
   label='CH4',
   reactive=True,
   structure=adjacencyList(
       """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
"""),
)

#----------
# Reaction systems
#1
surfaceReactor(
    temperature=(400,'K'),
    initialPressure=(15.0, 'bar'),
    initialGasMoleFractions={
        "CO": 0.2,
        "CO2": 0.2,
        "H2": 0.2,
        "N2": 0.2,
        "H2O": 0.2,  
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO2":0.15,},
    terminationTime=(1., 's'),
    terminationRateRatio=1e-3,
)
#2
surfaceReactor(
    temperature=(500,'K'),
    initialPressure=(15.0, 'bar'),
    initialGasMoleFractions={
        "CO": 0.2,
        "CO2": 0.2,
        "H2": 0.2,
        "N2": 0.2,
        "H2O": 0.2,  
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO2":0.15,},
    terminationTime=(1., 's'),
    terminationRateRatio=1e-3,
)
#3
surfaceReactor(
    temperature=(600,'K'),
    initialPressure=(15.0, 'bar'),
    initialGasMoleFractions={
        "CO": 0.2,
        "CO2": 0.2,
        "H2": 0.2,
        "N2": 0.2,
        "H2O": 0.2,  
    },
    initialSurfaceCoverages={
        "X": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    terminationConversion = { "CO2":0.15,},
    terminationTime=(10., 's'),
    terminationRateRatio=1e-3,
)

simulator(
    atol=1e-18,
    rtol=1e-12,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=500000,
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
    generateSeedEachIteration=True
)

generatedSpeciesConstraints(
    allowed=['input species','reaction libraries'],
    maximumOxygenAtoms=2,
    maximumCarbonAtoms=2,
    maximumSurfaceSites=2,
)

forbidden(
    label='H2Vdw',
    structure=adjacencyList(
        """
        1 H u0 p0 c0 {2,S}
        2 H u0 p0 c0 {1,S}
        3 X u0 p0 c0
        """
    )
)
forbidden(
    label="CO_bidentate_R",
        structure=adjacencyListGroup(
        """
        1 X u0 p0 c0 {2,S}
        2 O u0 p2 c0 {1,S} {3,S}
        3 C u0 p0 c0 {2,S} {4,[S,D,T]} {5,[S,D,T]}
        4 X u0 p0 c0 {3,[S,D,T]}
        5 R u0 c0 {3,[S,D,T]}
        """
    )
)
forbidden(
    label='ch4_vdw',
    structure=adjacencyList(
        """
        1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
        2 H u0 p0 c0 {1,S}
        3 H u0 p0 c0 {1,S}
        4 H u0 p0 c0 {1,S}
        5 H u0 p0 c0 {1,S}
        6 X u0 p0 c0
        """
    )
)
forbidden(
    label='cox2',
    structure=adjacencyList(
        """
        1 O u0 p2 c0 {2,S} {3,S}
        2 C u0 p0 c0 {1,S} {4,T}
        3 X u0 p0 c0 {1,S}
        4 X u0 p0 c0 {2,T}
        """
    )
)


