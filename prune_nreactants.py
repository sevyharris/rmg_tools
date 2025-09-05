import sys
from rmgpy.chemkin import load_chemkin_file, save_chemkin_file

def prune_reactions(input_chemkin, input_species_dict, output_chemkin, max_species=3):
    # Load chemkin file
    species_list, reaction_list = load_chemkin_file(input_chemkin, input_species_dict)
    
    # Filter reactions
    pruned_reactions = [
        rxn for rxn in reaction_list
        if len(rxn.reactants) <= max_species and len(rxn.products) <= max_species
    ]
    
    # Save pruned chemkin file
    save_chemkin_file(output_chemkin, species_list, pruned_reactions)

if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python prune_nreactants.py input.chemkin species_dictionary.txt output.chemkin [max_species]")
        sys.exit(1)
    input_chemkin = sys.argv[1]
    input_species_dict = sys.argv[2]
    output_chemkin = sys.argv[3]
    if len(sys.argv) == 5:
        max_species = int(sys.argv[4])
    else:
        max_species = 3
    prune_reactions(input_chemkin, input_species_dict, output_chemkin, max_species)