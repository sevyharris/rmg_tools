# script to create a dict of the cantera rxn index and the rmg index it corresponds to
import os
import sys
import pickle
import numpy as np
import concurrent.futures

import cantera as ct
import rmgpy.chemkin


MAX_WORKERS = 16
def same_reaction(rmg_rxn, ct_rxn):
    rmg_r = set([str(x) for x in rmg_rxn.reactants])
    rmg_p = set([str(x) for x in rmg_rxn.products])

    ct_r = set(ct_rxn.reactants.keys())
    ct_p = set(ct_rxn.products.keys())
    return rmg_r == ct_r and rmg_p == ct_p


# load the cantera and chemkin files: assumes they have the same base name and are in the same folder
chemkin = sys.argv[1]

working_dir = os.path.dirname(chemkin)
transport = os.path.join(working_dir, 'tran.dat')
species_dict = os.path.join(working_dir, 'species_dictionary.txt')
output_pickle_file = os.path.join(working_dir, 'ct2rmg_rxn.pickle')

# this takes a really long time to run, so try loading the file first
if os.path.exists(output_pickle_file):
    print(f'ck2rmg already exists {output_pickle_file}')
    exit(0)

species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(
    chemkin,
    dictionary_path=species_dict,
    transport_path=transport,
    use_chemkin_names=True
)
if '-gas' in chemkin:
    surface_path = os.path.join(working_dir, 'chem_annotated-surface.inp')
    species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(
    chemkin,
    dictionary_path=species_dict,
    transport_path=transport,
    surface_path=surface_path,
    use_chemkin_names=True
)

yaml_file = os.path.join(working_dir, 'chem_annotated.yaml')
if not os.path.exists(yaml_file):
    yaml_file = f'{chemkin[:-4]}.yaml'
if not os.path.exists(yaml_file):
    yaml_file = f'{chemkin[:-4]}.cti'
if not os.path.exists(yaml_file):
    print('no yaml file?')
    exit(1)

base_gas = ct.Solution(yaml_file)
surface = np.any([rmg_rxn.is_surface_reaction() for rmg_rxn in reaction_list])
if surface:
    base_gas, _ = ct.import_phases(yaml_file, ["gas", "surface1"])
    surf = ct.Interface(yaml_file, 'surface1')


def get_corresepondence(ct_index, ct_surf=False):

    if ct_surf:
        ct_surf_index = ct_index - base_gas.n_reactions
        start_guess = min(ct_index, len(reaction_list))

        # Assume from previous runs that the RMG index is always smaller than the Cantera index
        # so search backwards first from the current index
        for j in range(start_guess)[::-1]:
            if same_reaction(reaction_list[j], surf.reactions()[ct_surf_index]):
                return j

        # finish searching through the range
        for j in range(start_guess, len(reaction_list)):
            if same_reaction(reaction_list[j], surf.reactions()[ct_surf_index]):
                return j
    else:
        start_guess = min(ct_index, len(reaction_list))

        # Assume from previous runs that the RMG index is always smaller than the Cantera index
        # so search backwards first from the current index
        for j in range(start_guess)[::-1]:
            if same_reaction(reaction_list[j], base_gas.reactions()[ct_index]):
                return j

        # finish searching through the range
        for j in range(start_guess, len(reaction_list)):
            if same_reaction(reaction_list[j], base_gas.reactions()[ct_index]):
                return j

        # search through in regular order
        # for j in range(len(reaction_list)):
        #     if same_reaction(reaction_list[j], base_gas.reactions()[ct_index]):
        #         return j
    return -1


ct2rmg_rxn = {}




# do it in parallel
ct_indices = np.arange(0, base_gas.n_reactions)
is_surface_list = [False] * base_gas.n_reactions
if surface:
    ct_indices = np.arange(0, base_gas.n_reactions + surf.n_reactions)
    is_surface_list = [False] * base_gas.n_reactions + [True] * surf.n_reactions

with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
    for ct_index, rmg_index in zip(ct_indices, executor.map(
        get_corresepondence,
        ct_indices,
        is_surface_list,
    )):
        ct2rmg_rxn[ct_index] = rmg_index

# # do it in serial
# for i in range(len(base_gas.reactions())):
#     ct2rmg_rxn[i] = get_corresepondence(i)

# save the result
with open('ct2rmg_rxn.pickle', 'wb') as handle:
    pickle.dump(ct2rmg_rxn, handle, protocol=pickle.HIGHEST_PROTOCOL)
