import os
import sys
import copy
import rmgpy.data.rmg

# rmg_family_identifier.py

def identify_rmg_family(reaction):
    """
    Dummy function to determine RMG reaction family based on reactants and products.
    In a real scenario, this would use RMG's reaction family rules.
    """
    # load possible families from RMG database
    database = rmgpy.data.rmg.RMGDatabase()
    database.load(
        path = rmgpy.settings['database.directory'],
        thermo_libraries = ['primaryThermoLibrary',], # include  'surfaceThermoPt111' if doing all families
        transport_libraries = [],
        reaction_libraries = [],
        seed_mechanisms = [],
        kinetics_families = 'default',
        kinetics_depositories = [],
        depository = False,
    )
    for family in database.kinetics.families:
        if not database.kinetics.families[family].auto_generated:
            database.kinetics.families[family].add_rules_from_training(thermo_database=database.thermo)
            database.kinetics.families[family].fill_rules_by_averaging_up(verbose=True)

    family_matches = []
    labeled_rxns = {}
    for family in database.kinetics.families:
        try:
            new_reaction = copy.deepcopy(reaction)
            database.kinetics.families[family].add_atom_labels_for_reaction(new_reaction)
            template_labels = database.kinetics.families[family].get_reaction_template_labels(new_reaction)
            template = database.kinetics.families[family].retrieve_template(template_labels)
            kinetics = database.kinetics.families[family].get_kinetics_for_template(template, degeneracy=new_reaction.degeneracy)[0]
            family_matches.append(family)
            labeled_rxns[family] = new_reaction

        except (rmgpy.exceptions.ActionError, IndexError):
            pass
    return family_matches, labeled_rxns

if __name__ == "__main__":
    # Example usage
    # Replace these with actual reactant/product lists or SMILES strings

    # take in a reaction index from command line
    if len(sys.argv) > 1:
        reaction_index = int(sys.argv[1])
        sys.path.append(os.environ['DATABASE_DIR'])
        import database_fun

        reaction = database_fun.index2reaction(reaction_index)

    else:
        reaction = rmgpy.reaction.Reaction()
        reaction.reactants = [
            rmgpy.species.Species().from_smiles("C=C"),
            rmgpy.species.Species().from_smiles("[H]"),
        ]
        reaction.products = [
            rmgpy.species.Species().from_smiles("CH3"),
            rmgpy.species.Species().from_smiles("H2"),
        ]
    print("Reaction:", reaction)

    families, labeled_rxns = identify_rmg_family(reaction)
    print("Reaction families:", ", ".join(families))
    print("Labeled reactions:")
    for family in families:
        print('Reactants:')
        for sp in labeled_rxns[family].reactants:
            print(sp.molecule[0].get_all_labeled_atoms())
        print()
        print('Products:')
        for sp in labeled_rxns[family].products:
            print(sp.molecule[0].get_all_labeled_atoms())