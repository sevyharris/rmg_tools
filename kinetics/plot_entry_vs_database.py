import os
import rmgpy.data.kinetics
import rmgpy.data.rmg
import numpy as np
import sys
import copy
from ascii_plot import ASCIIPlot

def plot_kinetics(rxns, labels=None):
    """Function for plotting reaction kinetics
    Takes in a list of RMG reactions (rmgpy.reaction.Reaction) or a single reaction
    """
    markers = ['•', '*', '+', 'x', 'o', '◆', '■', '▲']
    if type(rxns) != list:
        rxns = [rxns]

    # Create ASCII plot
    plot = ASCIIPlot(width=100, height=30)
    plot.set_title('Reaction Kinetics')
    plot.set_xlabel('1000 / T (K^-1)')
    plot.set_ylabel('log10(k)')
    
    T = np.linspace(300, 3000, 1001)
    for m, rxn in enumerate(rxns):
        k = np.zeros(len(T))
        for i in range(0, len(T)):
            k[i] = rxn.get_rate_coefficient(T[i], 101325)
        
        # Prepare data for ASCII plot
        x_data = (1000.0 / T).tolist()
        y_data = np.log10(k).tolist()
        
        # Add plot with marker and label
        marker = markers[m % len(markers)]
        label = labels[m] if labels and m < len(labels) else f'Reaction {m+1}'
        plot.plot(x_data, y_data, marker=marker, label=label)
    
    # Display the plot
    plot.display()



# plot the kinetics entry and the corresponding reaction in the RMG database to compare
"""Usage:
python plot_entry.py <path_to_reactions.py> [entry_index]
"""
kinetics_lib_file = sys.argv[1]  # path to reactions.py
lib_dir = os.path.dirname(os.path.abspath(kinetics_lib_file))
lib_name = os.path.basename(lib_dir)
lib_container = os.path.dirname(lib_dir)

ark_kinetics_database = rmgpy.data.kinetics.KineticsDatabase()
ark_kinetics_database.load_libraries(lib_container)
print(f'{len(ark_kinetics_database.libraries[lib_name].entries)} entries loaded')

if len(sys.argv) == 3:
    index = int(sys.argv[2])
else:
    index = 0
rxn = ark_kinetics_database.libraries[lib_name].entries[index].item
rxn.kinetics = ark_kinetics_database.libraries[lib_name].entries[index].data

# now grab the kinetics from the database for the same reaction and plot them together
print(f'Loading RMG database and filling rules...')
database = rmgpy.data.rmg.RMGDatabase()
database.load(  # for now just load the 4 families we can compute
    path = rmgpy.settings['database.directory'],
    thermo_libraries = ['primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'DFT_QCI_thermo'],
    reaction_libraries = [],
    kinetics_families = ['H_Abstraction', 'R_Addition_MultipleBond', 'intra_H_migration', 'Disproportionation'],
)
for family in database.kinetics.families:
    if not database.kinetics.families[family].auto_generated:
        database.kinetics.families[family].add_rules_from_training(thermo_database=database.thermo)
        database.kinetics.families[family].fill_rules_by_averaging_up(verbose=True)

# get the corresponding family
rmg_rxn = rmgpy.reaction.Reaction(reactants=rxn.reactants, products=rxn.products)
for sp in rmg_rxn.reactants + rmg_rxn.products:
    sp.thermo = database.thermo.get_thermo_data(sp)


def reactions_in_same_direction(reaction1, reaction2):
    assert reaction1.is_isomorphic(reaction2), 'Reactions are not even isomorphic'
    if len(reaction1.reactants) != len(reaction2.reactants):
        return False
    reactants_to_match = [r for r in reaction1.reactants]
    counter = 0
    while reactants_to_match:
        for i in range(len(reaction2.reactants)):
            if reaction2.reactants[i].is_isomorphic(reactants_to_match[0]):
                reactants_to_match.remove(reactants_to_match[0])
                break
        else:
            return False

def get_reverse_reaction(reaction):
    assert reaction.kinetics is not None
    rev_reaction = copy.deepcopy(reaction)
    tmp_reactants = rev_reaction.reactants
    rev_reaction.reactants = rev_reaction.products
    rev_reaction.products = tmp_reactants

    if isinstance(reaction.kinetics, rmgpy.kinetics.arrhenius.ArrheniusEP):
        dHrxn = reaction.get_enthalpy_of_reaction(298)
        reaction.kinetics = reaction.kinetics.to_arrhenius(dHrxn)
    rev_reaction.kinetics = reaction.generate_reverse_rate_coefficient()
    return rev_reaction

for family in database.kinetics.families:
    try:
        database.kinetics.families[family].add_atom_labels_for_reaction(rmg_rxn)
        template_labels = database.kinetics.families[family].get_reaction_template_labels(rmg_rxn)
        template = database.kinetics.families[family].retrieve_template(template_labels)
        kinetics = database.kinetics.families[family].get_kinetics_for_template(template, degeneracy=rmg_rxn.degeneracy)[0]
        rmg_rxn.kinetics = kinetics
        break
    except rmgpy.exceptions.ActionError:
        continue
    else:
        print(f'Found matching template in family {family}: {template}')
        print(f'Kinetics from RMG database: {kinetics}')
        break

# make sure the reactions are in the same direction for plotting
rmg_label = f'RMG Database - {family}'
if not reactions_in_same_direction(rmg_rxn, rxn):
    rmg_rxn = get_reverse_reaction(rmg_rxn)
    rmg_label += ' (reversed)'

print(f'Plotting kinetics for reaction: {rxn}')
plot_kinetics([rxn, rmg_rxn], labels=[f'Entry {index}', rmg_label])
