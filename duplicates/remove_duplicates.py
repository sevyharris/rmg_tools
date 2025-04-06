# script to remove duplicates
import os
import sys
import rmgpy.chemkin


# load the mechanism
chemkin_file = sys.argv[1]
sp_dict = os.path.join(os.path.dirname(chemkin_file), 'species_dictionary.txt')
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, sp_dict, check_duplicates=False)


# count the bad duplicate pairs
bad_duplicate_pairs = []

for i in range(len(reaction_list)):
    for j in range(i):
        if reaction_list[i].is_isomorphic(reaction_list[j]):
            if not (reaction_list[i].duplicate and reaction_list[j].duplicate):
                bad_duplicate_pairs.append([i, j])


    if len(reaction_list[i].reactants) == 1 and len(reaction_list[i].products) == 1:
        # check for case of A (+M) <=> B (+M) not matched to 2A (+M) <=> 2B (+M) 
        
        other_reaction = rmgpy.reaction.Reaction()
        other_reaction.reactants = [reaction_list[i].reactants[0], reaction_list[i].reactants[0]]
        other_reaction.products = [reaction_list[i].products[0], reaction_list[i].products[0]]
        
        for j in range(len(reaction_list)):
            if j == i:
                continue
            if reaction_list[j].is_isomorphic(other_reaction):
                bad_duplicate_pairs.append([i, j])