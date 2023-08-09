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
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(
    chemkin,
    dictionary_path=species_dict,
    transport_path=transport,
    use_chemkin_names=True
)

if os.path.exists(f'{chemkin[:-4]}.yaml'):
    base_gas = ct.Solution(f'{chemkin[:-4]}.yaml')
elif os.path.exists(f'{chemkin[:-4]}.cti'):    
    base_gas = ct.Solution(f'{chemkin[:-4]}.cti')
else:
    print('no Cantera file?')
    exit(1)

output_pickle_file = os.path.join(working_dir, 'ct2rmg_rxn.pickle')


# this takes a really long time to run, so try loading the file first
if os.path.exists(output_pickle_file):
    print(f'ck2rmg already exists {output_pickle_file}')
    exit(0)


def get_corresepondence(ct_index):

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
ct_indices = np.arange(0, len(base_gas.reactions()))
with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
    for ct_index, rmg_index in zip(ct_indices, executor.map(
        get_corresepondence,
        ct_indices
    )):
        ct2rmg_rxn[ct_index] = rmg_index

# # do it in serial
# for i in range(len(base_gas.reactions())):
#     ct2rmg_rxn[i] = get_corresepondence(i)

# save the result
with open('ct2rmg_rxn.pickle', 'wb') as handle:
    pickle.dump(ct2rmg_rxn, handle, protocol=pickle.HIGHEST_PROTOCOL)
