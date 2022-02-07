import numpy as np
import adsorbate_thermo
import matplotlib.pyplot as plt


molecular_weight = 18.02
frequencies = [49.5, 68.6, 73.6, 102.0, 437.6, 452.9, 1596.3, 3675.6, 3787.0]


# Define the system parameters (previously read in from a file, soon to be read in again)
a = 3.912  # lattice constant for Pt
unit_cell_area = 9.0 * np.sqrt(3) / 4.0 * a ** 2.0
site_area = unit_cell_area / 9.0

my_calc = adsorbate_thermo.AdsorbateThermoCalc(molecular_weight, frequencies)
Q_trans, S_trans, Cp_trans, dH_trans = my_calc.get_translation_thermo()
# print(Q_trans)

Q_vib, S_vib, dH_vib, Cv_vib = my_calc.get_vibrational_thermo()
# print(Q_vib)

i = -2
print(my_calc.temperatures[i], Q_trans[i])

# check results with notebook so far
