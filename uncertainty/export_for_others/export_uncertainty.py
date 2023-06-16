# script to export mechanism uncertainty data to npy files

# Instructions:
# 1. change the database settings to match the input file
# 2. specify the mechanism path
# 3. make sure RMG-database is on the same commit as the one used to run RMG
# 4. make sure RMG-Py is on the correlated_uncertainty branch


import os
import numpy as np
# import rmgpy.tools.uncertainty_gao
import uncertainty_gao

# # Load the RMG Model
mech_dir = '/work/westgroup/harris.se/autoscience/reaction_calculator/models/base_rmg24'

chemkin = os.path.join(mech_dir, 'chem_annotated.inp')
species_dict = os.path.join(mech_dir, 'species_dictionary.txt')

# uncertainty = rmgpy.tools.uncertainty_gao.Uncertainty(output_directory=os.path.join(mech_dir, 'rmg_uncertainty'))
uncertainty = uncertainty_gao.Uncertainty(output_directory=os.path.join(mech_dir, 'rmg_uncertainty'))
uncertainty.load_model(chemkin, species_dict)

# --------------- CAUTION!!! Databases here must match the ones used to generate the mechanism
# note - this cell stalls out on Discovery
# these should be identical for the BASE RMG model whether 24 hour or one week
# The current git HEAD for RMG-database is:  # 24 hour
#     b'3e2cab5aa87a58c2a01f389c23be5b65cc649af8'
#     b'Tue Jul 19 13:58:30 2022 -0400'
# The current git HEAD for RMG-database is:  # 1 week
#     b'99d0678330ba22eedfdcae31eac50da33d1d5015'
#     b'Wed Aug 17 10:28:23 2022 -0400'

thermo_libs = [
    'BurkeH2O2',
    'primaryThermoLibrary',
    'FFCM1(-)',
    'CurranPentane',
    'Klippenstein_Glarborg2016',
    'thermo_DFT_CCSDTF12_BAC',
    'DFT_QCI_thermo',
    'CBS_QB3_1dHR',
]

kinetic_libs = [
    'FFCM1(-)',
    'CurranPentane',
    'combustion_core/version5',
    'Klippenstein_Glarborg2016',
    'BurkeH2O2inArHe',
    'BurkeH2O2inN2',
]

uncertainty.load_database(
    thermo_libraries=thermo_libs,
    kinetics_families='default',
    reaction_libraries=kinetic_libs,
    kinetics_depositories=['training'],
)


# Get the different kinetic and thermo sources
uncertainty.extract_sources_from_model()
uncertainty.assign_parameter_uncertainties()

np.save(os.path.join(mech_dir, 'gao_reaction_uncertainty.npy'), uncertainty.kinetic_input_uncertainties)
np.save(os.path.join(mech_dir, 'gao_species_uncertainty.npy'), uncertainty.thermo_input_uncertainties)
