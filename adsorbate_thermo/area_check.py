from ase.build import fcc111
import numpy as np

# a = 3.912
a = 4.0
surface = fcc111("Pt", (3, 3, 3), a=a, periodic=True)
print(surface)

a1 = surface.cell[0]
a2 = surface.cell[1]
print(a1)
print(a2)
print(np.linalg.norm(a2))

area = np.linalg.norm(np.cross(a1, a2))
print(area)

a = 3.912  # lattice constant for Pt
unit_cell_area = 9.0 * np.sqrt(3) / 4.0 * a ** 2.0
site_area = unit_cell_area / 9.0
