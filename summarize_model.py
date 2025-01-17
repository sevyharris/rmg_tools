# script to take in a chemkin or cantera model and print a short summary
import os
import sys
import rmgpy.chemkin
import cantera as ct


# read the file from argument
mech_file = sys.argv[1]

if mech_file.endswith('.inp'):
    # read in the chemkin file
    sp_dict = os.path.join(os.path.dirname(mech_file), 'species_dictionary.txt')
    species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(mech_file, sp_dict)
    print(f'Loaded Chemkin File: {mech_file}')
    print(f'Model contains {len(species_list)} species and {len(reaction_list)} reactions')
elif mech_file.endswith('.cti') or mech_file.endswitch('.yaml'):
    # read in the cantera file
    gas = ct.Solution(mech_file)
    print(f'Loaded Cantera File: {mech_file}')
    print(f'Model contains {gas.n_species} species and {gas.n_reactions} reactions')
else:
    raise ValueError('File must be a chemkin or cantera file')
