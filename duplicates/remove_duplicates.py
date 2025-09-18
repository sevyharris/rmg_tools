# script to remove duplicates -- not sure why it doesn't work
import os
import sys
import rmgpy.chemkin


# load the mechanism
chemkin_file = sys.argv[1]
sp_dict = os.path.join(os.path.dirname(chemkin_file), 'species_dictionary.txt')
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, sp_dict, check_duplicates=False)

verbose = False
if len(sys.argv) > 2:
    if sys.argv[2] in ['verbose=True', '--verbose', '--verbose=True']:
        verbose = True

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
            if reaction_list[j].is_isomorphic(other_reaction) and not (reaction_list[i].duplicate and reaction_list[j].duplicate):
                reaction_list[i].duplicate = True
                reaction_list[j].duplicate = True
                # bad_duplicate_pairs.append([i, j])


print(f'Found {len(bad_duplicate_pairs)} bad duplicate pairs')
print(bad_duplicate_pairs)

remove_list = []  # indices of reactions to remove
remove_list_items = []  # Reaction objects of the reactions to remove
repaired_list = []
for bad_duplicate_pair in bad_duplicate_pairs:
    i = bad_duplicate_pair[0]
    j = bad_duplicate_pair[1]

    # one PDEP and one library -- this is okay
    if isinstance(reaction_list[i], rmgpy.rmg.pdep.PDepReaction) and \
            isinstance(reaction_list[j], rmgpy.data.kinetics.library.LibraryReaction):
        # mark them as duplicates
        reaction_list[i].duplicate = True
        reaction_list[j].duplicate = True
        repaired_list.append([i, j])
    elif not isinstance(reaction_list[i], rmgpy.data.kinetics.library.LibraryReaction) and \
            isinstance(reaction_list[j], rmgpy.rmg.pdep.PDepReaction):
        # mark them as duplicates
        reaction_list[i].duplicate = True
        reaction_list[j].duplicate = True
        repaired_list.append([i, j])

    # one PDEP and one library -- this is okay
    if isinstance(reaction_list[i].kinetics, rmgpy.kinetics.chebyshev.Chebyshev) and \
            isinstance(reaction_list[j], rmgpy.data.kinetics.library.LibraryReaction):
        # mark them as duplicates
        reaction_list[i].duplicate = True
        reaction_list[j].duplicate = True
        repaired_list.append([i, j])
    elif not isinstance(reaction_list[i], rmgpy.data.kinetics.library.LibraryReaction) and \
            isinstance(reaction_list[j].kinetics, rmgpy.kinetics.chebyshev.Chebyshev):
        # mark them as duplicates
        reaction_list[i].duplicate = True
        reaction_list[j].duplicate = True
        repaired_list.append([i, j])

    # one library and one family -- this is bad. remove
    elif isinstance(reaction_list[i], rmgpy.data.kinetics.library.LibraryReaction) and \
            not isinstance(reaction_list[j], rmgpy.data.kinetics.library.LibraryReaction):
        remove_list.append(j)
        remove_list_items.append(reaction_list[j])
    elif not isinstance(reaction_list[i], rmgpy.data.kinetics.library.LibraryReaction) and \
            isinstance(reaction_list[j], rmgpy.data.kinetics.library.LibraryReaction):
        remove_list.append(i)
        remove_list_items.append(reaction_list[i])

    # one PDEP and one family -- this is okay
    elif isinstance(reaction_list[i], rmgpy.rmg.pdep.PDepReaction) and \
            not isinstance(reaction_list[j], rmgpy.rmg.pdep.PDepReaction):
        # mark them as duplicates
        reaction_list[i].duplicate = True
        reaction_list[j].duplicate = True
        repaired_list.append([i, j])
    elif not isinstance(reaction_list[i], rmgpy.rmg.pdep.PDepReaction) and \
            isinstance(reaction_list[j], rmgpy.rmg.pdep.PDepReaction):
        # mark them as duplicates
        reaction_list[i].duplicate = True
        reaction_list[j].duplicate = True
        repaired_list.append([i, j])

    # two PDEPs -- this is probably wrong
    elif isinstance(reaction_list[i], rmgpy.rmg.pdep.PDepReaction) and \
            isinstance(reaction_list[j], rmgpy.rmg.pdep.PDepReaction):
        # remove the one with more reactants
        if len(reaction_list[i].reactants) > reaction_list[j].reactants:
            remove_list.append(i)
            remove_list_items.append(reaction_list[i])
        elif len(reaction_list[i].reactants) < reaction_list[j].reactants:
            remove_list.append(j)
            remove_list_items.append(reaction_list[j])
        else:  # remove j by default because it was added second
            remove_list.append(j)
            remove_list_items.append(reaction_list[j])

    else:
        print('Unknown kind of bad pair -- not implemented')
        print(type(reaction_list[i]))
        print(type(reaction_list[j]))
        print()

# actually remove the items
for i, item in enumerate(remove_list_items):
    reaction_list.remove(item)
    if verbose:
        print(f'Removing reaction {remove_list[i]}: {item}')

if verbose:
    for pair in repaired_list:
        print(f'Repaired {pair}')

print(f'Removed {len(remove_list_items)}')
print(f'Repaired {len(repaired_list)}')

# save the resulting chemkin file
out_chemkin = chemkin_file.replace('.inp', '_fixed.inp')
rmgpy.chemkin.save_chemkin_file(out_chemkin, species_list, reaction_list)
