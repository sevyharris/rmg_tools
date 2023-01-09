# script to take in a chemkin or cantera model and print a short summary

import sys
import rmgpy.chemkin
import cantera as ct


# read the file from argument
mech_file = sys.argv[1]

if mech_file.endswith('.inp'):
    # read in the chemkin file
    species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(mech_file)
    print(f'Loaded Chemkin File: {mech_file}')
    print(f'Model contains {len(species_list)} species and {len(reaction_list)} reactions')
elif mech_file.endswith('.cti'):
    # read in the cantera file
    gas = ct.Solution(mech_file)
    print(f'Loaded Cantera File: {mech_file}')
    print(f'Model contains {gas.n_species} species and {gas.n_reactions} reactions')
else:
    raise ValueError('File must be a chemkin or cantera file')
