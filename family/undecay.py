# script to figure out the original reaction before rmg did the decay recipe in rmgpy.rmg.decay

import sys

import rmgpy.data.rmg
import rmgpy.species
import rmgpy.rmg.decay

sys.path.append('/work/westgroup/harris.se/autoscience/reaction_calculator/dft')
import autotst_wrapper

reaction_index = 4752
reaction_index = 4732
rmg_reaction = autotst_wrapper.database_fun.index2reaction(reaction_index)

database = rmgpy.data.rmg.RMGDatabase()
decay_families = ['intra_H_migration']  # Add to this if you find more examples

database.load(
    path = rmgpy.settings['database.directory'],
    thermo_libraries = ['primaryThermoLibrary'],
    transport_libraries = [],
    reaction_libraries = [],
    seed_mechanisms = [],
    kinetics_families = decay_families,
    kinetics_depositories = ['training'],
    #frequenciesLibraries = self.statmechLibraries,
    depository = False,
)
for family in database.kinetics.families:
    if not database.kinetics.families[family].auto_generated:
        database.kinetics.families[family].add_rules_from_training(thermo_database=database.thermo)
        database.kinetics.families[family].fill_rules_by_averaging_up(verbose=True)

mol_reaction = rmgpy.reaction.Reaction()
mol_reaction.reactants = [sp.molecule[0] for sp in rmg_reaction.reactants]
mol_reaction.products = [sp.molecule[0] for sp in rmg_reaction.products]

possible_reactions = []
for family in decay_families:
    possible_reactions += database.kinetics.families[family].generate_reactions([sp.molecule[0] for sp in rmg_reaction.reactants])

decayed_matches = []
for i in range(len(possible_reactions)):
    
    new_reaction = rmgpy.reaction.Reaction()
    new_reaction.reactants = rmg_reaction.reactants
    
    decayed_products = []
    for p in possible_reactions[i].products:
        decayed_products += rmgpy.rmg.decay.decay_species(rmgpy.species.Species(molecule=[p]))
    
    new_reaction.products = decayed_products
    if new_reaction.is_isomorphic(rmg_reaction):
        decayed_matches.append(possible_reactions[i])

for d in decayed_matches:
    print(d)

for d in decayed_matches:
    query_reaction = rmgpy.reaction.Reaction()
    query_reaction.reactants = rmg_reaction.reactants
    query_reaction.products = [rmgpy.species.Species(molecule=[m]) for m in d.products]
    
    print(d, autotst_wrapper.database_fun.get_unique_reaction_index(query_reaction))
