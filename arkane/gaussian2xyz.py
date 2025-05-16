# script to convert a Gaussian .log file to an ase .xyz gile

import ase
import sys
import arkane.ess.gaussian


def get_xyz(atoms, fmt='%22.15f'):  # taken from ase.io.xyz's write_xyz function
    # returns the positions of the ase.atoms.Atoms object in xyz format as a string
    xyz_str = ""
    for s, (x, y, z) in zip(atoms.symbols, atoms.positions):
        xyz_str += ('%-2s %s %s %s\n' % (s, fmt % x, fmt % y, fmt % z))
    return xyz_str


input_path = sys.argv[1]
gaussian_logfile = arkane.ess.gaussian.GaussianLog(input_path, check_for_errors=False)

coord, number, mass = gaussian_logfile.load_geometry()
atoms = ase.Atoms(positions=coord, symbols=number)
print(get_xyz(atoms))
