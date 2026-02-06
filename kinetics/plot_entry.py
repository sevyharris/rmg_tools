import os
import rmgpy.data.kinetics
import numpy as np
import sys
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
print(f'Plotting kinetics for reaction: {rxn}')
plot_kinetics(rxn, labels=[f'Entry {index}'])


