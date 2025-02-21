import ase.build
import numpy as np
import rmgpy.constants

################ Pt(111) #####################
print('------------------------- Pt(111) -----------------')
a = 3.912
surface = ase.build.fcc111("Pt", (3, 3, 3), a=a, periodic=True)

a1 = surface.cell[0]
a2 = surface.cell[1]
unit_cell_area = np.linalg.norm(np.cross(a1, a2))
print(f'Unit cell area:\t{unit_cell_area:0.5f}\tA^2')

site_area = unit_cell_area / 9.0
print(f'Site area:\t{site_area:0.5f}\t\tA^2/site')
print(f'Site density:\t{1.0 / site_area:0.5f}\t\tsites/A^2')

SDEN = 1.0 / site_area * (1e20) / rmgpy.constants.Na
print(f'Site density:\t{SDEN:0.5e}\tmols/m^2')
SDEN = 1.0 / site_area * (1e16) / rmgpy.constants.Na
print(f'Site density:\t{SDEN:0.5e}\tmols/cm^2')

#a = 3.912  # lattice constant for Pt
#unit_cell_area = 9.0 * np.sqrt(3) / 4.0 * a ** 2.0
#site_area = unit_cell_area / 9.0

################## Fe(110) ###################
print()
print('------------------------- Fe(110) -----------------')
a = 2.85  # A
surface = ase.build.bcc110("Fe", (3, 3, 3), a=a, periodic=True)

a1 = surface.cell[0]
a2 = surface.cell[1]
unit_cell_area = np.linalg.norm(np.cross(a1, a2))
print(f'Unit cell area:\t{unit_cell_area:0.5f}\tA^2')

site_area = unit_cell_area / 9.0
print(f'Site area:\t{site_area:0.5f}\t\tA^2/site')
print(f'Site density:\t{1.0 / site_area:0.5f}\t\tsites/A^2')

SDEN = 1.0 / site_area * (1e20) / rmgpy.constants.Na
print(f'Site density:\t{SDEN:0.5e}\tmols/m^2')
SDEN = 1.0 / site_area * (1e16) / rmgpy.constants.Na
print(f'Site density:\t{SDEN:0.5e}\tmols/cm^2')

