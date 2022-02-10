# ethane ranked list of uncertain/sensitive parameters
import os
from rmgpy.tools.uncertainty import Uncertainty, process_local_results
from rmgpy.tools.canteramodel import get_rmg_species_from_user_species
from rmgpy.species import Species


# Load the model
model_dir = '/home/moon/rmg/my_examples/ethane/'
chemkin_file = os.path.join(model_dir, 'chemkin', 'chem_annotated.inp')
dict_file = os.path.join(model_dir, 'chemkin', 'species_dictionary.txt')

# Initialize the Uncertainty class instance and load the model
uncertainty = Uncertainty(output_directory=os.path.join(model_dir, 'chemkin'))
uncertainty.load_model(chemkin_file, dict_file)


# sensitivity
# Map the species to the objects within the Uncertainty class
ethane = Species().from_smiles('CC')
C2H4 = Species().from_smiles('C=C')
mapping = get_rmg_species_from_user_species([ethane, C2H4], uncertainty.species_list)

# Define the reaction conditions
initial_mole_fractions = {mapping[ethane]: 1.0}
T = (1300, 'K')
P = (1, 'bar')
termination_time = (0.5, 'ms')
sensitive_species = [mapping[ethane], mapping[C2H4]]

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

print('sensitivity complete')

uncertainty.load_database(
    thermo_libraries=['primaryThermoLibrary'],
    kinetics_families='default',
    reaction_libraries=[],
)

uncertainty.extract_sources_from_model()
uncertainty.assign_parameter_uncertainties()

result = uncertainty.local_analysis(sensitive_species, correlated=False, number=5, fileformat='.png')
print(process_local_results(result, sensitive_species, number=5)[1])
