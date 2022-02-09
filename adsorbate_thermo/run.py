import numpy as np
import adsorbate_thermo
import rmgpy.species
import arkane.output
import rmgpy.data.thermo


molecular_weight = 18.02
frequencies = [49.5, 68.6, 73.6, 102.0, 437.6, 452.9, 1596.3, 3675.6, 3787.0]
composition = {'H': 2, 'O': 1, 'C': 0, 'N': 0, 'Pt': 1}
heat_of_formation_0K = -259.05  # kJ/mol

# Define the system parameters (previously read in from a file, soon to be read in again)
# could also get from cross product of unit cell vectors in ase
a = 3.912  # lattice constant for Pt
unit_cell_area = 9.0 * np.sqrt(3) / 4.0 * a ** 2.0
site_area = unit_cell_area / 9.0

my_calc = adsorbate_thermo.AdsorbateThermoCalc(molecular_weight, frequencies, composition, heat_of_formation_0K, twoD_gas=True)
h2o = rmgpy.species.Species(smiles='O')
h2o.thermo = my_calc.get_thermo()

species_list = [h2o]
path = '/home/moon/rmg/rmg_tools/adsorbate_thermo/'
name = 'sevy2'
lib_long_desc = 'sevy''s second attempt at exporting adsorbate thermo'

arkane.output.save_thermo_lib(species_list, path, name, lib_long_desc)
