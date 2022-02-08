import numpy as np
import adsorbate_thermo
import matplotlib.pyplot as plt
import rmgpy.species

# name = 'H2O_ads'
# DFT_binding_energy = [-1.89E-01, 'eV']
# heat_of_formation_0K = [-259.05, 'kJ/mol']
# composition = {'H':2, 'C':0, 'N':0, 'O':1, 'Pt':1}
# sites = 1
# adsorbate_mass = [18.02, 'amu']
# linear_scaling_binding_atom = [-3.5860E+00, 'eV']
# linear_scaling_gamma(X) = [0]
# linear_scaling_psi = [-0.189318, 'eV']
# frequencies = [49.5, 68.6, 73.6, 102.0, 437.6, 452.9, 1596.3, 3675.6, 3787.0, 'cm-1']
# ~

molecular_weight = 18.02
frequencies = [49.5, 68.6, 73.6, 102.0, 437.6, 452.9, 1596.3, 3675.6, 3787.0]
composition = {'H': 2, 'O': 1, 'C': 0, 'N': 0}
heat_of_formation_0K = -259.05  # kJ/mol

# Define the system parameters (previously read in from a file, soon to be read in again)
a = 3.912  # lattice constant for Pt
unit_cell_area = 9.0 * np.sqrt(3) / 4.0 * a ** 2.0
site_area = unit_cell_area / 9.0

my_calc = adsorbate_thermo.AdsorbateThermoCalc(molecular_weight, frequencies, composition, heat_of_formation_0K, twoD_gas=True)
Q_trans, S_trans, Cp_trans, dH_trans = my_calc.get_translation_thermo()
# print(Q_trans)

Q_vib, S_vib, dH_vib, Cv_vib = my_calc.get_vibrational_thermo()
# print(Q_vib)

my_calc.get_thermo()
# i = 0
# # print(my_calc.temperatures[i], Q_vib[i])
# # print(my_calc.temperatures[i], S_vib[i])
# print(my_calc.temperatures[i], dH_vib[i])
# # print(my_calc.temperatures[i], Cv_vib[i])
# print(dH_vib.shape)
# # check results with notebook so far
