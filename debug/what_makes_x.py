# script to see what makes a given species
import rmgpy.chemkin
import rmgpy.species

import rmgpy.tools.fluxdiagram

import numpy as np
import os
import sys
sys.path.append(os.environ['DATABASE_DIR'])
import database_fun

# Define input species labels (edit these as needed)
smiles1 = 'CC1OCC1'   # 449 <--- not found in the mechanism

# Paths to your Chemkin and species dictionary files
chemkin_path = '/home/moon/autoscience/compare_to_aramco/branching/twenty_separate/rmg_min_3/chem_annotated.inp'
chemkin_path =  '/home/moon/autoscience/compare_to_aramco/branching/concatenate/chem_G1p5_lnk1p3.inp'
chemkin_path =  '/home/moon/autoscience/compare_to_aramco/branching/twenty_separate/rmg_min_4_bigger/chem_annotated.inp'
chemkin_path =  '/home/moon/autoscience/aramco/pruned_aramco.inp'
chemkin_path = '/home/moon/autoscience/fuels/butane_20240501/chem_annotated.inp'

species_dict_path = os.path.join(os.path.dirname(chemkin_path), 'species_dictionary.txt')


# Load species and reactions from Chemkin file
# species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_path, species_dict_path)


# --------------------------------------------------------
# simulate the reaction to figure out the highest flux options
rmg_input_file = '/home/moon/autoscience/fuels/simplified_input.py'  # for conditions T = 830, P = 10 atm
t_end = 0.1  # the ignition time for Aramco

print('Loading RMG job 1...')
rmg_job1 = rmgpy.tools.fluxdiagram.load_rmg_job(
    rmg_input_file,
    chemkin_path,
    species_dict_path,
    generate_images=False,
    check_duplicates=True
)

rmg_job1.reaction_systems[0].termination[0].time.value_si = t_end  # override the termination time

print('Conducting simulation')
times1, concentrations1, reaction_rates1 = rmgpy.tools.fluxdiagram.simulate(
    rmg_job1.reaction_model,
    rmg_job1.reaction_systems[0],

)
print('Simulation complete')  # termination is set by the rmg_input_file, but you could probably override it here too
species_list1 = rmg_job1.reaction_model.core.species[:]
reaction_list1 = rmg_job1.reaction_model.core.reactions[:]

ref_sp1 = rmgpy.species.Species(smiles=smiles1)

ref_index1 = database_fun.get_unique_species_index(ref_sp1)
print(f"Looking for reactions involving species indices: {ref_index1} ({smiles1})")

# Find reactions where one input species goes to another input species
matching_reaction_indices = []
for j, rxn in enumerate(reaction_list1):
    if np.any([ref_sp1.is_isomorphic(rxn.reactants[i]) for i in range(len(rxn.reactants))]) or \
       np.any([ref_sp1.is_isomorphic(rxn.products[i]) for i in range(len(rxn.products))]):
        matching_reaction_indices.append(j)

# for each reaction in the matching reactions, find its cumulative net flux over the simulation time
reaction_fluxes = np.zeros(len(matching_reaction_indices), float)
for j in range(len(matching_reaction_indices)):
    reaction_fluxes[j] = np.sum(reaction_rates1[:, matching_reaction_indices[j]])


# Sort reactions by flux
indices = np.arange(len(matching_reaction_indices))
sorted_order = [x for _, x in sorted(zip(reaction_fluxes, indices))][::-1]  # descending order

# list the top 20 reactions
for i in range(20):
    if i < len(sorted_order):
        percent_flux = reaction_fluxes[sorted_order[i]] / np.sum(reaction_fluxes) * 100
        try:
            db_index = database_fun.get_unique_reaction_index(reaction_list1[matching_reaction_indices[sorted_order[i]]])
        except IndexError:
            db_index = -1
        print(f"{db_index}\t{percent_flux:.2f}\t{reaction_list1[matching_reaction_indices[sorted_order[i]]]}")


# # Print matching reactions
# print(f"Reactions where one input species goes to another input species:")
# if not matching_reaction_indices:
#     print("No matching reactions found.")
# else:
#     for j in matching_reaction_indices:
#         print(reaction_list[j])
