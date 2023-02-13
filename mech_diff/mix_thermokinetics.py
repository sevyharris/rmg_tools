# A module with functions to mix two reaction kinetics or species thermos

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import cantera as ct

import rmgpy.chemkin
import rmgpy.thermo.thermodata
import rmgpy.quantity


def mix_kinetics(reaction0, reaction1, w):
    # function to mix reaction 0 and reaction 1 according to weight (a slider from 0-1)
    # 0 is all reaction 0 and 1 is all reaction 1
    # this method fits the kinetics to the a set of points a fraction w of the way between reactions 0 and 1

    # note, you'll get better results if you just change the parameters a little at a time...

    # Use the RMG thermo object rmgpy.thermo.nasa.NASA

    kinetics0 = reaction0.kinetics
    kinetics1 = reaction1.kinetics

    if w == 0.0:
        return kinetics0
    elif w == 1.0:
        return kinetics1

    # print(type(kinetics0))
    assert type(kinetics0) == rmgpy.kinetics.arrhenius.Arrhenius
    assert type(kinetics1) == rmgpy.kinetics.arrhenius.Arrhenius
    assert reaction0.kinetics.A.units == reaction1.kinetics.A.units
    assert reaction0.kinetics.n.units == reaction1.kinetics.n.units
    assert reaction0.kinetics.Ea.units == reaction1.kinetics.Ea.units

    N = 1001
    As = np.logspace(
        np.log10(reaction0.kinetics.A.value),
        np.log10(reaction1.kinetics.A.value),
        N
    )

    ns = np.linspace(
        reaction0.kinetics.n.value,
        reaction1.kinetics.n.value,
        N
    )

    Eas = np.linspace(
        reaction0.kinetics.Ea.value,
        reaction1.kinetics.Ea.value,
        N
    )

    i = int(w * (N - 1))
    mixed_rxn = rmgpy.kinetics.arrhenius.Arrhenius(
        A=(As[i], reaction0.kinetics.A.units),
        n=(ns[i], reaction0.kinetics.n.units),
        Ea=(Eas[i], reaction0.kinetics.Ea.units)
    )

    return mixed_rxn


def mix_kinetics_retired(reaction0, reaction1, w):
    # function to mix reaction 0 and reaction 1 according to weight (a slider from 0-1)
    # 0 is all reaction 0 and 1 is all reaction 1

    # Use the RMG thermo object rmgpy.thermo.nasa.NASA

    kinetics0 = reaction0.kinetics
    kinetics1 = reaction1.kinetics

    if w == 0.0:
        return kinetics0
    elif w == 1.0:
        return kinetics1

    # print(type(kinetics0))
    assert type(kinetics0) == rmgpy.kinetics.arrhenius.Arrhenius
    assert type(kinetics1) == rmgpy.kinetics.arrhenius.Arrhenius

    try:
        Tmin = np.maximum(kinetics0.Tmin.value, kinetics1.Tmin.value)
    except AttributeError:
        Tmin = 300
    try:
        Tmax = np.minimum(kinetics0.Tmax.value, kinetics1.Tmax.value)
    except AttributeError:
        Tmax = 3000

    N_total = 1001
    Tdata = np.linspace(Tmin, Tmax, N_total)
    kdata = np.zeros(N_total)
    lnkdata = np.zeros(N_total)
    P = 101325.0
    for i, T in enumerate(Tdata):
        lnkdata[i] = w * np.log(reaction0.get_rate_coefficient(T, P)) + (1.0 - w) * np.log(reaction1.get_rate_coefficient(T, P))
        kdata[i] = np.exp(lnkdata[i])

    # get the units
    # kunits = rmgpy.kinetics.model.get_rate_coefficient_units_from_reaction_order(len(reaction1.reactants))
    kunits = reaction1.kinetics.A.units

    mixed_kinetics = rmgpy.kinetics.arrhenius.Arrhenius()
    mixed_kinetics.fit_to_data(Tdata, kdata, kunits)

    return mixed_kinetics


def mix_thermo(species0, species1, w):
    # function to mix species 0 and species 1 according to weight (a slider from 0-1)
    # 0 is all species 0 and 1 is all species 1

    # Use the RMG thermo object rmgpy.thermo.nasa.NASA

    thermo0 = species0.thermo
    thermo1 = species1.thermo

    if w == 0.0:
        return thermo0
    elif w == 1.0:
        return thermo0

    Tmin = np.maximum(thermo0.Tmin.value, thermo1.Tmin.value)
    Tmax = np.minimum(thermo0.Tmax.value, thermo1.Tmax.value)

    N_total = 1001
    Tdata = np.linspace(Tmin, Tmax, N_total)
    Cpdata = np.zeros(N_total)
    for i, T in enumerate(Tdata):
        Cpdata[i] = w * thermo0.get_heat_capacity(T) + (1.0 - w) * thermo1.get_heat_capacity(T)

    # generate all the Tdata and Cpdata
    # put 298 in first position - I know, this is erasing one of the Tdata's but I don't care
    # ThermoData is requiring that this be greater than 298 for some reason
    Tdata[0] = 300.0
    Cpdata[0] = w * thermo0.get_heat_capacity(Tdata[0]) + (1.0 - w) * thermo1.get_heat_capacity(Tdata[0])

    Tdata = rmgpy.quantity.Quantity(Tdata, 'K')
    Cpdata = rmgpy.quantity.Quantity(Cpdata, 'J/(mol*K)'),

    # Always use H298 and S298 from the first thermo given
    H298 = w * thermo0.get_enthalpy(298.0) + (1.0 - w) * thermo1.get_enthalpy(298.0)
    H298 = rmgpy.quantity.Quantity(H298, 'J/mol')
    S298 = w * thermo0.get_entropy(298.0) + (1.0 - w) * thermo1.get_entropy(298.0)
    S298 = rmgpy.quantity.Quantity(S298, 'J/(mol*K)')

    Cp0 = w * species0.calculate_cp0() + (1.0 - w) * species1.calculate_cp0()
    Cp0 = rmgpy.quantity.Quantity(Cp0, 'J/(mol*K)')
    CpInf = w * species0.calculate_cpinf() + (1.0 - w) * species1.calculate_cpinf()
    CpInf = rmgpy.quantity.Quantity(CpInf, 'J/(mol*K)')

    my_data = rmgpy.thermo.thermodata.ThermoData(
        Tdata=Tdata,
        Cpdata=Cpdata,
        H298=H298,
        S298=S298,
        Cp0=Cp0,
        CpInf=CpInf,
    )
    Tint = 1000.0
    if Tint > Tmax or Tint < Tmin:
        Tint = (Tmax + Tmin) / 2.0
    nasa = my_data.to_nasa(Tmin, Tmax, Tint)
    return nasa


def plot_kinetics(kinetics, labels=None):
    plt.xlabel('1000 / T (K^-1)')
    plt.ylabel('ln(k)')

    T = np.linspace(300, 3000, 1001)
    for kinetic in kinetics:
        k = np.zeros(len(T))
        for i in range(0, len(T)):
            k[i] = kinetic.get_rate_coefficient(T[i])
        plt.plot(1000.0 / T, np.log(k))

    plt.legend(labels)
    plt.show()


def plot_thermos(thermos, labels=None):
    """
    Function to plot a list of rmgpy.thermo.nasa.NASA (RMG reaction.thermo) objects
    """
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
    for thermo in thermos:
        H = np.zeros(len(T))
        S = np.zeros(len(T))
        Cp = np.zeros(len(T))
        for i in range(0, len(T)):
            H[i] = thermo.get_enthalpy(T[i]) / 1000.0
            S[i] = thermo.get_entropy(T[i]) / 1000.0
            Cp[i] = thermo.get_heat_capacity(T[i]) / 1000.0
        ax[0].plot(T, H)
        ax[1].plot(T, S)
        ax[2].plot(T, Cp)
    ax[0].legend(labels)
    ax[1].legend(labels)
    ax[2].legend(labels)
    plt.subplots_adjust(wspace=0.25)
    plt.show()


def test_thermo_mix():
    print('Loading ethane thermo example')

    matplotlib.use('TkAgg')

    # load an rmg/chemkin model just for reaction/species objects to play with
    chemkin_path = '/home/moon/rmg/my_examples/ethane/chemkin/chem_annotated.inp'
    dictionary_path = '/home/moon/rmg/my_examples/ethane/chemkin/species_dictionary.txt'
    transport_path = '/home/moon/rmg/my_examples/ethane/chemkin/tran.dat'
    species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_path, dictionary_path=dictionary_path, transport_path=transport_path)

    sp4 = species_list[4]
    sp5 = species_list[5]

    thermo4p5 = mix_thermo(sp4, sp5, w=0.5)
    plot_thermos(
        [sp4.thermo, sp5.thermo, thermo4p5], labels=[str(sp4), str(sp5), 'mix']
    )


def test_kinetics_mix():
    print('Loading ethane kinetics example')

    matplotlib.use('TkAgg')

    # load an rmg/chemkin model just for reaction/species objects to play with
    chemkin_path = '/home/moon/rmg/my_examples/ethane/chemkin/chem_annotated.inp'
    dictionary_path = '/home/moon/rmg/my_examples/ethane/chemkin/species_dictionary.txt'
    transport_path = '/home/moon/rmg/my_examples/ethane/chemkin/tran.dat'
    species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_path, dictionary_path=dictionary_path, transport_path=transport_path)

    rxn4 = reaction_list[4]
    rxn5 = reaction_list[5]

    rxn4p5 = mix_kinetics(rxn4, rxn5, w=0.5)
    plot_kinetics(
        [rxn4.kinetics, rxn5.kinetics, rxn4p5], labels=[str(rxn4), str(rxn5), 'mix']
    )


if __name__ == '__main__':
    test_thermo_mix()
    test_kinetics_mix()
