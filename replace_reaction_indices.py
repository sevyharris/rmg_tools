import argparse
import os
import numpy as np
import rmgpy.chemkin
import sys
import copy
sys.path.append(os.environ['DATABASE_DIR'])
import database_fun


def parse_comma_separated_integers(arg_string):
    try:
        return [int(item.strip()) for item in arg_string.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid comma-separated integers format.")


def get_i_thing(thing, thing_list, pdep_specific=False):
    if not pdep_specific:
        for i in range(len(thing_list)):
            if thing.is_isomorphic(thing_list[i]):
                return i
        return -1
    
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            if (isinstance(thing, rmgpy.rmg.pdep.PDepReaction) and isinstance(thing_list[i], rmgpy.rmg.pdep.PDepReaction)) or \
                    (not isinstance(thing, rmgpy.rmg.pdep.PDepReaction) and not isinstance(thing_list[i], rmgpy.rmg.pdep.PDepReaction)):
                return i
    return -1


def get_all_i_things(thing, thing_list):
    indices = []
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            indices.append(i)
    return indices

def main():
    parser = argparse.ArgumentParser(description="Replace reactions in a Chemkin file with those from a reference Chemkin file.")
    parser.add_argument("chemkin_file", help="Path to the input Chemkin file")
    parser.add_argument("reference_chemkin_file", help="Path to the reference Chemkin file")
    parser.add_argument("-s", "--species", type=parse_comma_separated_integers, help="Indices (comma separated) of species to replace")
    parser.add_argument("-r", "--reactions", type=parse_comma_separated_integers, help="Indices (comma separated) of reactions to replace")
    parser.add_argument("-o", "--output", default="chemkin_replaced.inp", help="Output Chemkin file name")
    args = parser.parse_args()

    species_dict = os.path.join(os.path.dirname(args.chemkin_file), 'species_dictionary.txt')
    reference_species_dict = os.path.join(os.path.dirname(args.reference_chemkin_file), 'species_dictionary.txt')

    # Load original chemkin file
    species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(args.chemkin_file, species_dict)
    # Load reference chemkin file
    ref_species_list, ref_reaction_list = rmgpy.chemkin.load_chemkin_file(args.reference_chemkin_file, reference_species_dict)

    # Replace specified reactions
    sp_added = False
    for idx in args.reactions:
        ref_reaction = database_fun.index2reaction(idx)
        i_rxn_refs = get_all_i_things(ref_reaction, ref_reaction_list)
        i_rxn_news = sorted(get_all_i_things(ref_reaction, reaction_list))[::-1]  # reverse sort so we can delete without messing up indices

        if len(i_rxn_news) <= 0:
            print(f"Reaction with index {idx} not found in input chemkin file. Will just add it from reference chemkin.")
            

        # assert len(i_rxn_news) > 0, f"Reaction with index {idx} not found in input chemkin file."


        # in either case, we will delete all the reactions from the original and replace with all the reactions from the reference
        if len(i_rxn_refs) <= 0:
            print(f"Reaction with index {idx} not found in reference chemkin file. Deleting all references from input chemkin.")
            for i_rxn in i_rxn_news:
                reaction_list.pop(i_rxn)
            continue

        # delete all the reactions from the original
        for i_rxn in i_rxn_news:
            reaction_list.pop(i_rxn)

        # now add in all the reactions from the reference
        for i_rxn_ref in i_rxn_refs:
            # make a copy of the reference reaction but using the species from the original species list
            # new_reaction = copy.deepcopy(ref_reaction_list[i_rxn_ref])
            new_reaction = ref_reaction_list[i_rxn_ref]


            # add the species
            for i in range(len(new_reaction.reactants)):
                j = get_i_thing(new_reaction.reactants[i], species_list)
                if j < 0:
                    print(f"Species {new_reaction.reactants[i]} not found in input chemkin file. Adding it.")
                    species_list.append(new_reaction.reactants[i])
                    j = len(species_list) - 1
                    sp_added = True
                new_reaction.reactants[i] = species_list[j]
            for i in range(len(new_reaction.products)):
                j = get_i_thing(new_reaction.products[i], species_list)
                if j < 0:
                    print(f"Species {new_reaction.products[i]} not found in input chemkin file. Adding it.")
                    species_list.append(new_reaction.products[i])
                    j = len(species_list) - 1
                    sp_added = True
                new_reaction.products[i] = species_list[j]

            new_reaction.generate_pairs()
            reaction_list.append(new_reaction)
            print(f"Added reaction {new_reaction} from reference chemkin file.")
            print(new_reaction.kinetics)
            print()
            print()

    reactions_to_remove = []
    for idx in args.species:
        ref_species = database_fun.index2species(idx)
        i_sp_ref = get_i_thing(ref_species, ref_species_list)
        i_sp_new = get_i_thing(ref_species, species_list)

        assert i_sp_new >= 0, f"Species with index {idx} not found in input chemkin file."

        if i_sp_ref < 0:
            print(f"Species with index {idx} not found in reference chemkin file. Deleting it from input chemkin.")
            species_list.pop(i_sp_new)

            # todo we also need to remove any reactions that reference this species
            for i_rxn, rxn in enumerate(reaction_list):
                if np.any([ref_species.is_isomorphic(rxn.reactants[i]) for i in range(len(rxn.reactants))]) or np.any([ref_species.is_isomorphic(rxn.products[i]) for i in range(len(rxn.products))]):
                    reactions_to_remove.append(i_rxn)
            continue

        # replace the thermo data
        species_list[i_sp_new].thermo = ref_species_list[i_sp_ref].thermo

    reactions_to_remove = sorted(set(reactions_to_remove))[::-1]  # reverse sort so we can delete without messing up indices
    print(f"Removing {len(reactions_to_remove)} reactions that reference deleted species.")
    # Remove reactions that reference deleted species
    for i_rxn in reactions_to_remove:
        reaction_list.pop(i_rxn)

    # Save new chemkin file
    rmgpy.chemkin.save_chemkin_file(args.output, species_list, reaction_list)
    if sp_added:
        rmgpy.chemkin.save_species_dictionary(species_dict, species_list)

if __name__ == "__main__":
    main()


# sample usage:
# python ~/rmg/rmg_tools/replace_reaction_indices.py /home/moon/autoscience/compare_to_aramco/branching/twenty_separate/rmg_min_3/chem_annotated.inp /home/moon/autoscience/compare_to_aramco/branching/concatenate/chem_G1p5_lnk1p3.inp -s 4,5 -r 2470,4721 -o chemkin_replaced.inp
