# change reactor settings to single point so the flux diagram generator doesn't crash

# Data sources
database(
    thermoLibraries = ['BurkeH2O2', 'CurranPentane','FFCM1(-)','primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo', 'CBS_QB3_1dHR'],
    reactionLibraries = ['CurranPentane','FFCM1(-)','combustion_core/version5'],
    seedMechanisms = ['BurkeH2O2inN2','BurkeH2O2inArHe', 'C2H4+O_Klipp2017'],#added BurkeH2O2inAr and BurkeH2O2inN2 for bath gases
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)
# List of species 
species(
    label = 'nheptane',
    reactive = True,
    structure = SMILES("CCCCCCC"),
)

species(
    label = 'O2',
    reactive = True,
    structure = SMILES("[O][O]"),
)

species(
    label = 'N2',
    reactive = False,
    structure = SMILES("N#N"),
)

species(
    label = 'Ar',
    reactive = False,
    structure = SMILES("[Ar]"),
)


# We want CH, C2H and O to be forced into the core,
# and to have known (reproducible) species numbers
# because we use them when post-processing ignition
# delay simulations.

species(
    label = 'CH',
    reactive = True,
    structure = SMILES("[CH]"),
)

species(
    label = 'C2H',
    reactive = True,
    structure = SMILES("[C]#C"),
)

species(
    label = 'O',
    reactive = True,
    structure = SMILES("[O]"),
)


#Species that typically show up in JSR reactors based on Dagaut 2014
species(
    label = 'CO',
    reactive = True,
    structure = SMILES("[C-]#[O+]"),
)

species(
    label = 'CO2',
    reactive = True,
    structure = SMILES("O=C=O"),
)

species(
    label = 'H2O',
    reactive = True,
    structure = SMILES("O"),
)

species(
    label = 'CH2O',
    reactive = True,
    structure = SMILES("C=O"),
)

species(
    label = 'CH4',
    reactive = True,
    structure = SMILES("C"),
)

species(
    label = 'C2H4',
    reactive = True,
    structure = SMILES("C=C"),
)

species(
    label = 'C3H6',
    reactive = True,
    structure = SMILES("C=CC"),
)

species(
    label = 'H2',
    reactive = True,
    structure = SMILES("[H][H]"),
)


# Reaction systems 1 for the high temperature range with Argon as the bath gas 
simpleReactor(
    temperature=(550,'K'),
    # temperature=[(550,'K'),(1500,'K')],#max reactor temperature must always be less than max pdep temperature
    pressure=(1.0,'bar'),
    # pressure=[(1.0,'bar'),(40.0,'bar')],
    nSims=5,
    initialMoleFractions={
        # "nheptane": 0.003747, # phi = 1
        "nheptane": 0.003747, # range of 0.5 < phi < 2
        # "nheptane": [0.003747/2, 0.003747*2], # range of 0.5 < phi < 2
        "Ar": 0.955040,
        "O2": 0.041213,
        },
    terminationTime=(1.0, 's'),
    terminationRateRatio=0.01,
)

# # Reaction systems 2 for the high temperature range with Nitrogen as the bath gas 
# simpleReactor(
#     temperature=[(550,'K'),(1500,'K')],#max reactor temperature must always be less than max pdep temperature
#     pressure=[(1.0,'bar'),(40.0,'bar')],
#     nSims=5,
#     initialMoleFractions={
#         # "nheptane": 0.003747, # phi = 1
#         "nheptane": [0.003747/2, 0.003747*2], # range of 0.5 < phi < 2
#         "N2": 0.955040,
#         "O2": 0.041213,
#         },
#     terminationTime=(1.0, 's'),
#     terminationRateRatio=0.01,
# )

# # Reaction system 3 for the low temperature range with Argon as the bath gas
# simpleReactor(
#     temperature=[(550,'K'),(800,'K')],
#     pressure=[(1.0,'bar'),(40.0,'bar')],
#     nSims=17,
#     initialMoleFractions={
#         # "nheptane": 0.003747, # phi = 1
#         "nheptane": [ 0.003747/2, 0.003747*2], # range of 0.5 < phi < 2
#         "Ar": 0.955040,
#         "O2": 0.041213,
#         },
#     terminationTime=(1.0, 's'),
#     terminationRateRatio=0.01,
# )


# # Reaction system 4 for the low temperature range with Nitrogen as the bath gas
# simpleReactor(
#     temperature=[(550,'K'),(800,'K')],
#     pressure=[(1.0,'bar'),(40.0,'bar')],
#     nSims=17,
#     initialMoleFractions={
#         # "nheptane": 0.003747, # phi = 1
#         "nheptane": [ 0.003747/2, 0.003747*2], # range of 0.5 < phi < 2
#         "N2": 0.955040,
#         "O2": 0.041213,
#         },
#     terminationTime=(1.0, 's'),
#     terminationRateRatio=0.01,
# )

# Introducing staging into the model and simulator blocks

# first stage with a higher tolerance
simulator(
    atol=1e-16,
    rtol=1e-8,
)

# # second stage with a lower tolerance 
# simulator(
#     atol=1e-16,
#     rtol=1e-8,
# )

# first stage with a lower maxNumSpecies=300
model(
    toleranceKeepInEdge=0.0055,
    toleranceMoveToCore=0.1,
    filterReactions = True,
    filterThreshold = 5e8, # Default is 5e8
    maximumEdgeSpecies=100000,
    maxNumObjsPerIter=1,# number of objects (species, reactions, PDepNetworks) that can be taken from one simulation
    maxNumSpecies = 300,
)


# # first stage with a lower maxNumSpecies=500

# model(
#     toleranceKeepInEdge=0.0055,
#     toleranceMoveToCore=0.1,
#     filterReactions = True,
#     filterThreshold = 5e8, # Default is 5e8
#     maximumEdgeSpecies=100000,
#     maxNumObjsPerIter=3,
#     maxNumSpecies = 500,
# )

pressureDependence(
        method='modified strong collision',
        maximumGrainSize=(0.5,'kcal/mol'),
        minimumNumberOfGrains=250,
        temperatures=(300,2000,'K',8),#due to thermo library max pdep=2500K cannot go to 3000K
        pressures=(0.01,100,'bar',5),
        interpolation=('Chebyshev', 6, 4),
        maximumAtoms=16,
)

options(
    units='si',
    saveRestartPeriod=None,
    generateOutputHTML=True,
    generatePlots=False,
    saveEdgeSpecies=True,
    saveSimulationProfiles=True,
 )
# Introduced the species constraint for the mechanism. Set the maximumRadicalElectrons to 2 for N2
generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries'],
    maximumCarbonAtoms=10,
    maximumRadicalElectrons=2,
)

