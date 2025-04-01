# useful functions for plotting

import numpy as np
import matplotlib.pyplot as plt


def plot_kinetics(rxns, labels=None):
    """Function for plotting reaction kinetics
    Takes in a list of RMG reactions (rmgpy.reaction.Reaction) or a single reaction
    """
    plt.xlabel('1000 / T (K^-1)')
    plt.ylabel('log10(k)')

    if type(rxns) != list:
        rxns = [rxns]

    T = np.linspace(300, 3000, 1001)
    for rxn in rxns:
        k = np.zeros(len(T))
        for i in range(0, len(T)):
            k[i] = rxn.get_rate_coefficient(T[i], 101325)
        plt.plot(1000.0 / T, np.log10(k))

    if labels:
        plt.legend(labels)
    plt.show()


def plot_thermos(thermos, labels=None):
    if type(thermos) != list:
        thermos = [thermos]
    if labels is None:
        labels = ['' for t in thermos]
    linestyles = ['solid', 'dashed', 'dotted']
    fig, ax = plt.subplots(1, 3)
    fig.set_size_inches(12, 3)
    fig.tight_layout()
    ax[0].set_xlabel('Temperature (K)')
    ax[0].set_ylabel('H (kJ / mol)')
    ax[0].set_title('Enthalpy vs. Temperature')
    ax[1].set_xlabel('Temperature (K)')
    ax[1].set_ylabel('S (kJ / mol K)')
    ax[1].set_title('Entropy vs. Temperature')
    ax[2].set_xlabel('Temperature (K)')
    ax[2].set_ylabel('Cp (kJ / mol K)')
    ax[2].set_title('Heat Capacity vs. Temperature')
    T = np.linspace(300, 3000, 1001)
    for m, thermo in enumerate(thermos):
        if 'cantera' in str(type(thermo)).lower() and 'species' in str(type(thermo)).lower():
            thermo = thermo.thermo
        H = np.zeros(len(T))
        S = np.zeros(len(T))
        Cp = np.zeros(len(T))
        if 'rmgpy' in str(type(thermo)).lower():
            for i in range(0, len(T)):
                H[i] = thermo.get_enthalpy(T[i]) / 1000.0
                S[i] = thermo.get_entropy(T[i]) / 1000.0
                Cp[i] = thermo.get_heat_capacity(T[i]) / 1000.0
        else:  # cantera
            for i in range(0, len(T)):
                H[i] = thermo.h(T[i]) / 1e6  # J/kmol
                S[i] = thermo.s(T[i]) / 1e6
                Cp[i] = thermo.cp(T[i]) / 1e6  # J/mol K
        ax[0].plot(T, H, linestyle=linestyles[m % len(linestyles)])
        ax[1].plot(T, S, linestyle=linestyles[m % len(linestyles)])
        ax[2].plot(T, Cp, linestyle=linestyles[m % len(linestyles)])
    ax[0].legend(labels)
    ax[1].legend(labels)
    ax[2].legend(labels)
    ax[2].yaxis.get_major_formatter().set_useOffset(False)
    plt.subplots_adjust(wspace=0.25)
    plt.show()
