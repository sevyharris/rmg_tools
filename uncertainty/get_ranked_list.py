import os
import tabulate
import pandas as pd
from rmgpy.tools.uncertainty import Uncertainty, process_local_results
from rmgpy.tools.canteramodel import get_rmg_species_from_user_species
from rmgpy.species import Species


# Load the model
# Must use annotated chemkin file

model_dir = '/home/moon/rmg/rmg_tools/uncertainty/nheptane'


chemkin_file = os.path.join(model_dir, 'chem_annotated.inp')
dict_file = os.path.join(model_dir, 'species_dictionary.txt')

# Initialize the Uncertainty class instance and load the model
uncertainty = Uncertainty(output_directory=os.path.join(model_dir, 'chemkin'))
uncertainty.load_model(chemkin_file, dict_file)

# Map the species to the objects within the Uncertainty class
nheptane = Species().from_smiles('CCCCCCC')
CO2 = Species().from_smiles('O=C=O')
O2 = Species().from_smiles('[O][O]')
OH_rad = Species().from_smiles('[OH]')
Ne = Species().from_smiles('[Ne]')
mapping = get_rmg_species_from_user_species([nheptane, CO2, O2, OH_rad, Ne], uncertainty.species_list)


# TODO get model with He instead of Ne
# equivalence ratio phi, 11 O2 per n-heptane in stoichiometric feed
x_nheptane = 0.005
phi = 1.0
x_O2 = 11.0 * x_nheptane / phi
x_Ne = 1.0 - x_nheptane - x_O2
# print(x_nheptane)
# print(x_O2)
# print(x_Ne)

# Reaction conditions to match at least some Zhang(2016) experimental data
initial_mole_fractions = {
    mapping[nheptane]: x_nheptane,
    mapping[O2]: x_O2,
    mapping[Ne]: x_Ne,
}
T = (1000, 'K')
P_torr = 800.0
P_atm = P_torr / 760.0
P_Pa = P_atm * 101325.0
P_bar = P_Pa / 100000.0

P = (P_bar, 'bar')
termination_time = (0.5, 'ms')
sensitive_species = [mapping[nheptane], mapping[CO2], mapping[O2], mapping[OH_rad]]

# Perform the sensitivity analysis
uncertainty.sensitivity_analysis(
    initial_mole_fractions,
    sensitive_species,
    T,
    P,
    termination_time,
    number=5,
    fileformat='.png'
)

# NOTE: You must load the database with the same settings which were used to generate the model.
#       This includes any thermo or kinetics libraries which were used.
uncertainty.load_database(
    thermo_libraries=[
        'BurkeH2O2',
        'CurranPentane',
        'FFCM1(-)',
        'primaryThermoLibrary',
        'thermo_DFT_CCSDTF12_BAC',
        'DFT_QCI_thermo',
        'CBS_QB3_1dHR'
    ],
    kinetics_families='default',
    reaction_libraries=[
        'CurranPentane',
        'FFCM1(-)',
        'combustion_core/version5'
    ]
)

uncertainty.extract_sources_from_model()
uncertainty.assign_parameter_uncertainties()

result = uncertainty.local_analysis(sensitive_species, correlated=False, number=15, fileformat='.png')
print(process_local_results(result, sensitive_species, number=15)[1])

