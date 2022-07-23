import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import rmgpy.thermo
import cantera as ct


def ct2rmg_thermo(ct_thermo):
    rmg_thermo = rmgpy.thermo.nasa.NASA()
    rmg_thermo.Tmin = (ct_thermo.min_temp, 'K')
    rmg_thermo.Tmax = (ct_thermo.max_temp, 'K')
    midpoint = ct_thermo.coeffs[0]
    a_high = ct_thermo.coeffs[1:8]
    a_low = ct_thermo.coeffs[8:]

    nasa1 = rmgpy.thermo.nasa.NASAPolynomial()
    nasa2 = rmgpy.thermo.nasa.NASAPolynomial()

    nasa1.coeffs = a_low
    nasa1.Tmin = rmg_thermo.Tmin
    nasa1.Tmax = (midpoint, 'K')

    nasa2.coeffs = a_high
    nasa2.Tmin = (midpoint, 'K')
    nasa2.Tmax = rmg_thermo.Tmax

    rmg_thermo.polynomials = [nasa1, nasa2]

    return rmg_thermo


def test_ct2rmg_thermo():
    # load the gri30 mech for an example
    print("Running Cantera to RMG thermo")
    gas = ct.Solution('gri30.yaml')
    ct_thermo = gas.species()[0].thermo
    rmg_thermo = ct2rmg_thermo(ct_thermo)
    # print(ct_thermo)
    # print(rmg_thermo)
    # print(rmg_thermo.get_enthalpy(500))

    # plot the thermo results
    fig, ax = plt.subplots(1, 3)
    fig.set_size_inches(12, 3)
    fig.tight_layout()
    ax[0].set_xlabel('Temperature (K)')
    ax[0].set_ylabel('H (J / mol)')
    ax[0].set_title('Enthalpy vs. Temperature')
    ax[1].set_xlabel('Temperature (K)')
    ax[1].set_ylabel('S (J / mol K)')
    ax[1].set_title('Entropy vs. Temperature')
    ax[2].set_xlabel('Temperature (K)')
    ax[2].set_ylabel('Cp (J / mol K)')
    ax[2].set_title('Heat Capacity vs. Temperature')

    T = np.linspace(300, 3000, 1001)

    H_rmg = np.zeros(len(T))
    S_rmg = np.zeros(len(T))
    Cp_rmg = np.zeros(len(T))
    H_ct = np.zeros(len(T))
    S_ct = np.zeros(len(T))
    Cp_ct = np.zeros(len(T))
    for i in range(0, len(T)):
        H_rmg[i] = rmg_thermo.get_enthalpy(T[i])
        S_rmg[i] = rmg_thermo.get_entropy(T[i])
        Cp_rmg[i] = rmg_thermo.get_heat_capacity(T[i])
        H_ct[i] = ct_thermo.h(T[i]) / 1000.0
        S_ct[i] = ct_thermo.s(T[i]) / 1000.0
        Cp_ct[i] = ct_thermo.cp(T[i]) / 1000.0
    labels = ['RMG', 'Cantera']
    ax[0].plot(T, H_rmg, T, H_ct, ':')
    ax[1].plot(T, S_rmg, T, S_ct, ':')
    ax[2].plot(T, Cp_rmg, T, Cp_ct, ':')
    ax[0].legend(labels)
    ax[1].legend(labels)
    ax[2].legend(labels)
    plt.subplots_adjust(wspace=0.25)
    plt.show()


if __name__ == '__main__':
    test_ct2rmg_thermo()
