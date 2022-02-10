# Adapted from Notebook from Prof. C. Franklin Goldsmith and Katrín Blöndal at Brown University

import numpy as np
import rmgpy.constants


R = rmgpy.constants.R
kB = rmgpy.constants.kB
h = rmgpy.constants.h
c = rmgpy.constants.c
amu = rmgpy.constants.amu
Avogadro = rmgpy.constants.Na
pi = np.pi
invcm_to_invm = 1.0E2


T_switch = 1000.0


class AdsorbateThermoCalc:
    """
    Class for molecule stat thermo

    site area - assumes Pt(111) fcc facet
    """
    def __init__(
        self,
        molecular_weight,           # g/mol
        frequencies,                # TODO make a load function, pretty sure this needs to be sorted
        composition,                # TODO put this in load function
        heat_of_formation_0K,       # TODO put this in load function
        site_area=None,             # surface area per binding site m^2, default is Pt(111) facet
        site_occupation_number=1,   # number of sites occupied by adsorbate
        cutoff_frequency=100.0,     # cm^-1
        twoD_gas=False,
        temperatures=None,
    ):
        self.molecular_weight = molecular_weight
        self.frequencies = frequencies
        self.composition = composition
        self.heat_of_formation_0K = heat_of_formation_0K
        self.site_occupation_number = site_occupation_number
        if site_area is None:
            site_area = 62.10456e-20 / 9.0
        self.site_area = site_area
        self.cutoff_frequency = cutoff_frequency
        self.twoD_gas = twoD_gas
        if temperatures is None:
            # NOTE 298.15 must be first for the NASA polynomial routine to work!
            temperatures = [298.15]
            T_low = 300.0
            T_high = 2000.0
            dT = 10.0  # temperature increment
            temperatures = np.append(temperatures, np.arange(T_low, T_high + dT, dT))
        self.temperatures = temperatures

    def load_molecule_info(self):
        pass

    def get_Cp0(self):
        # TODO verify this
        # units = h * c / kB * invcm_to_invm

        # frequencies = self.frequencies  # TODO check sorted
        # if self.twoD_gas:  # skip the first two if we do 2D gas
        #     frequencies = frequencies[2:]

        # Cv_vib = 0
        # temp = 0.0
        # t = 0
        # for (n, nu) in enumerate(self.frequencies):
        #     # x goes to infinity
        #     x = 1e30
        #     Cv_vib[t] += x**2.0 * np.exp(-x) / (1.0 - np.exp(-x)) ** 2.0

        # Cv_vib *= R
        # return 0
        return rmgpy.quantity.Quantity(R, 'kJ/(mol*K)')

    def get_CpInf(self):
        # used wolfram alpha to get limits
        # TODO verify this
        frequencies = self.frequencies  # TODO check sorted
        if self.twoD_gas:  # skip the first two if we do 2D gas
            frequencies = frequencies[2:]

        return rmgpy.quantity.Quantity(R * len(frequencies), 'kJ/(mol*K)')

    def get_translation_thermo(self):
        # Just using the area is not a function of temperature equations
        Q_trans = ((2.0 * pi * self.molecular_weight * amu * kB * self.temperatures)
                   / np.float_power(h, 2.0)) * self.site_area * self.site_occupation_number
        S_trans = R * (2.0 + np.log(Q_trans))
        Cp_trans = R * np.ones(len(self.temperatures))  # NOTE: Cp = Cv
        dH_trans = R * 1.0 * self.temperatures
        return Q_trans, S_trans, Cp_trans, dH_trans

    def get_vibrational_thermo(self):
        units = h * c / kB * invcm_to_invm

        frequencies = self.frequencies  # TODO check sorted
        if self.twoD_gas:  # skip the first two if we do 2D gas
            frequencies = frequencies[2:]

        # TODO point to a LaTeX notebook explaining the equations AND the matrixization
        # x is a matrix (Temperatures * Frequencies in size)
        # x = nu * units / temp #cm^-1 * K cm / K = dimensionless
        x = np.matmul(np.matrix(units / self.temperatures).transpose(), np.matrix(frequencies))
        Q_vib = np.prod(1.0 / (1.0 - np.exp(-x)), 1)
        S_vib = np.sum(-np.log(1.0 - np.exp(-x)) + np.multiply(x, np.exp(-x)) / (1.0 - np.exp(-x)), 1) * R
        dH_vib = np.multiply(np.sum(np.multiply(x, np.exp(-x)) / (1.0 - np.exp(-x)), 1), np.matrix(self.temperatures).transpose()) * R
        Cv_vib = np.sum(np.multiply(np.float_power(x, 2.0), np.exp(-x)) / np.float_power((1.0 - np.exp(-x)), 2.0), 1) * R

        # rewrap the matrices as arrays
        Q_vib = np.squeeze(np.array(Q_vib))
        S_vib = np.squeeze(np.array(S_vib))
        dH_vib = np.squeeze(np.array(dH_vib))
        Cv_vib = np.squeeze(np.array(Cv_vib))
        return Q_vib, S_vib, dH_vib, Cv_vib

    def get_thermo(self):
        # call the subroutine for the vibrational partition function
        Q_trans, S_trans, Cp_trans, dH_trans = self.get_translation_thermo()
        Q_vib, S_vib, dH_vib, Cv_vib = self.get_vibrational_thermo()

        # now compute the correction to the heat of formation as you go from 0 to 298 K
        # TODO figure out where this came from
        h_correction = 4.234  # kJ/mol. enthalpy_H(298) - enthalpy_H(0)
        c_correction = 1.051  # kJ/mol. enthalpy_C(298) - enthalpy_C(0)
        n_correction = 4.335  # kJ/mol. enthalpy_N(298) - enthalpy_N(0)
        o_correction = 4.340  # kJ/mol. enthalpy_O(298) - enthalpy_O(0)
        heat_of_formation_correction = 0.0
        heat_of_formation_correction += self.composition['H'] * h_correction
        heat_of_formation_correction += self.composition['C'] * c_correction
        heat_of_formation_correction += self.composition['N'] * n_correction
        heat_of_formation_correction += self.composition['O'] * o_correction

        # note that the partition function is the production of the individual terms,
        # whereas the thermodynamic properties are additive
        Q = Q_trans * Q_vib
        S = S_trans + S_vib
        dH = dH_trans + dH_vib
        Cp = Cp_trans + Cv_vib  # see comments in each section regarding Cp vs Cv
        heat_of_formation_298K = self.heat_of_formation_0K + dH[0] - heat_of_formation_correction
        H = heat_of_formation_298K + dH - dH[0]
        self.S1 = S

        index298 = list(self.temperatures).index(298.15)

        import rmgpy.quantity
        Tdata = rmgpy.quantity.Quantity(self.temperatures, 'K')
        Cpdata = rmgpy.quantity.Quantity(Cp, 'kJ/(mol*K)'),
        H298 = rmgpy.quantity.Quantity(H[index298], 'kJ/mol'),
        S298 = rmgpy.quantity.Quantity(S[index298], 'kJ/(mol*K)'),

        Cp0 = self.get_Cp0()
        CpInf = self.get_CpInf()

        Tmin = np.min(self.temperatures)
        Tmax = np.max(self.temperatures)
        Tint = 1000.0
        import rmgpy.thermo.thermodata
        my_data = rmgpy.thermo.thermodata.ThermoData(
            Tdata=Tdata,
            Cpdata=Cpdata,
            H298=H298,
            S298=S298,
            Cp0=Cp0,
            CpInf=CpInf,
        )
        nasa = my_data.to_nasa(Tmin, Tmax, Tint)
        return nasa

    def get_thermo2(self):
        # katrin's notebook
        # call the subroutine for the vibrational partition function
        Q_trans, S_trans, Cp_trans, dH_trans = self.get_translation_thermo()
        Q_vib, S_vib, dH_vib, Cv_vib = self.get_vibrational_thermo()

        # now compute the correction to the heat of formation as you go from 0 to 298 K
        # TODO figure out where this came from
        h_correction = 4.234  # kJ/mol. enthalpy_H(298) - enthalpy_H(0)
        c_correction = 1.051  # kJ/mol. enthalpy_C(298) - enthalpy_C(0)
        n_correction = 4.335  # kJ/mol. enthalpy_N(298) - enthalpy_N(0)
        o_correction = 4.340  # kJ/mol. enthalpy_O(298) - enthalpy_O(0)
        heat_of_formation_correction = 0.0
        heat_of_formation_correction += self.composition['H'] * h_correction
        heat_of_formation_correction += self.composition['C'] * c_correction
        heat_of_formation_correction += self.composition['N'] * n_correction
        heat_of_formation_correction += self.composition['O'] * o_correction

        # note that the partition function is the production of the individual terms,
        # whereas the thermodynamic properties are additive
        Q = Q_trans * Q_vib
        S = S_trans + S_vib
        dH = dH_trans + dH_vib
        Cp = Cp_trans + Cv_vib  # see comments in each section regarding Cp vs Cv
        heat_of_formation_298K = self.heat_of_formation_0K + dH[0] - heat_of_formation_correction
        H = heat_of_formation_298K + dH - dH[0]
        self.S2 = S

        a_low, a_high = self.fit_NASA2(Cp, H, S)
        return a_low, a_high

    def fit_NASA2(self, Cp, H, S,):

        heat_capacity = Cp
        reference_enthalpy = H[0]
        reference_entropy = S[0]

        i_switch = -1
        for i in range(len(self.temperatures)):
            if self.temperatures[i] == T_switch:
                i_switch = i
        if i_switch == -1:
            print("We have a problem! Cannot find switching temperature")

        # start by creating the independent variable matrix for the low-temperature fit
        YT = np.array([np.ones(len(self.temperatures[:i_switch + 1])), self.temperatures[:i_switch + 1], self.temperatures[:i_switch + 1] ** 2.0, self.temperatures[:i_switch + 1] ** 3.0, self.temperatures[:i_switch + 1] ** 4.0], dtype=np.float64)  # this is transpose of our Y
        Y = YT.transpose()  # this is the desired Y

        b = heat_capacity[:i_switch + 1] / R
        a_low = np.linalg.lstsq(Y, b, rcond=None)[0]

        T_ref = 298.15
        # now determine the enthalpy coefficient for the low-T region
        subtract = a_low[0] + a_low[1] / 2.0 * T_ref + a_low[2] / 3.0 * T_ref**2.0 + a_low[3] / 4.0 * T_ref**3.0 + a_low[4] / 5.0 * T_ref**4.0
        a_low = np.append(a_low, reference_enthalpy / R - subtract * T_ref)
        # now determine the entropy coefficient for the low-T region
        subtract = a_low[0] * np.log(T_ref) + a_low[1] * T_ref + a_low[2] / 2.0 * T_ref**2.0 + a_low[3] / 3.0 * T_ref**3.0 + a_low[4] / 4.0 * T_ref**4.0
        a_low = np.append(a_low, reference_entropy / R - subtract)

        #
        # NOW SWITCH TO HIGH-TEMPERATURE REGIME!
        #
        T_ref = T_switch
        # compute the heat capacity, enthalpy, and entropy at the switching point
        Cp_switch = a_low[0] + a_low[1] * T_ref + a_low[2] * T_ref**2.0 + a_low[3] * T_ref**3.0 + a_low[4] * T_ref**4.0
        H_switch = a_low[0] * T_ref + a_low[1] / 2.0 * T_ref**2.0 + a_low[2] / 3.0 * T_ref**3.0 + a_low[3] / 4.0 * T_ref**4.0 + a_low[4] / 5.0 * T_ref**5.0 + a_low[5]
        S_switch = a_low[0] * np.log(T_ref) + a_low[1] * T_ref + a_low[2] / 2.0 * T_ref**2.0 + a_low[3] / 3.0 * T_ref**3.0 + a_low[4] / 4.0 * T_ref**4.0 + a_low[6]

        # now repeat the process for the high-temperature regime
        a_high = [0.0]
        YT = np.array([self.temperatures[i_switch:], self.temperatures[i_switch:]**2.0, self.temperatures[i_switch:]**3.0, self.temperatures[i_switch:]**4.0], dtype=np.float64)  # this is transpose of our Y
        Y = YT.transpose()  # this is the desired Y

        b = heat_capacity[i_switch:] / R - Cp_switch
        a_high = np.append(a_high, np.linalg.lstsq(Y, b, rcond=None)[0])
        a_high[0] = Cp_switch - (a_high[0] + a_high[1] * T_switch + a_high[2] * T_switch**2.0 + a_high[3] * T_switch**3.0 + a_high[4] * T_switch**4.0)

        a_high = np.append(a_high, H_switch - (a_high[0] + a_high[1] / 2.0 * T_ref + a_high[2] / 3.0 * T_ref**2.0 + a_high[3] / 4.0 * T_ref**3.0 + a_high[4] / 5.0 * T_ref**4.0) * T_ref)
        a_high = np.append(a_high, S_switch - (a_high[0] * np.log(T_ref) + a_high[1] * T_ref + a_high[2] / 2.0 * T_ref**2.0 + a_high[3] / 3.0 * T_ref**3.0 + a_high[4] / 4.0 * T_ref**4.0))

        # Check to see if there is a discontinuity
        if (1 == 0):
            print("\ncheck for discontinuities:")
            cp_low_Tswitch = a_low[0] + a_low[1] * T_switch + a_low[2] * T_switch**2.0 + a_low[3] * T_switch**3.0 + a_low[4] * T_switch**4.0
            cp_high_Tswitch = a_high[0] + a_high[1] * T_switch + a_high[2] * T_switch**2.0 + a_high[3] * T_switch**3.0 + a_high[4] * T_switch**4.0
            H_low_Tswitch = a_low[0] * T_switch + a_low[1] / 2.0 * T_switch**2.0 + a_low[2] / 3.0 * T_switch**3.0 + a_low[3] / 4.0 * T_switch**4.0 + a_low[4] / 5.0 * T_switch**5.0 + a_low[5]
            H_high_Tswitch = a_high[0] * T_switch + a_high[1] / 2.0 * T_switch**2.0 + a_high[2] / 3.0 * T_switch**3.0 + a_high[3] / 4.0 * T_switch**4.0 + a_high[4] / 5.0 * T_switch**5.0 + a_high[5]
            S_low_Tswitch = a_low[0] * np.log(T_switch) + a_low[1] * T_switch + a_low[2] / 2.0 * T_switch**2.0 + a_low[3] / 3.0 * T_switch**3.0 + a_low[4] / 4.0 * T_switch**4.0 + a_low[6]
            S_high_Tswitch = a_high[0] * np.log(T_switch) + a_high[1] * T_switch + a_high[2] / 2.0 * T_switch**2.0 + a_high[3] / 3.0 * T_switch**3.0 + a_high[4] / 4.0 * T_switch**4.0 + a_high[6]

            print("discontinuity at T_switch for Cp/R is %.4F" % (cp_low_Tswitch - cp_high_Tswitch))
            print("discontinuity at T_switch for H/R is %.4F" % (H_low_Tswitch - H_high_Tswitch))
            print("discontinuity at T_switch for S/R is %.4F" % (S_low_Tswitch - S_high_Tswitch))

        # line = '\n\t !cut and paste this value into the cti file!\n'
        line = '\tthermo = (\n'
        line += "\t\tNASA( [%.1F, %.1F], [%.8E, %.8E,\n \t\t %.8E, %.8E, %.8E,\n \t\t %.8E, %.8E]), \n" % (300.0, 1000.0, a_low[0], a_low[1], a_low[2], a_low[3], a_low[4], a_low[5], a_low[6])
        line += "\t\tNASA( [%.1F, %.1F], [%.8E, %.8E,\n \t\t %.8E, %.8E, %.8E,\n \t\t %.8E, %.8E]), \n" % (1000.0, max(self.temperatures), a_high[0], a_high[1], a_high[2], a_high[3], a_high[4], a_high[5], a_high[6])
        line += "\t\t ),\n"

        # molecule.thermo_lines = line
        return a_low, a_high
        # molecule.a_low = a_low
        # molecule.a_high = a_high

    def get_thermo_from_NASA2(self, a_low, a_high):
        # compute thermo properties from nasa polynomials

        i_switch = -1
        for i in range(len(self.temperatures)):
            if self.temperatures[i] == T_switch:
                i_switch = i

        cp_fit = np.zeros(len(self.temperatures))
        h_fit = np.zeros(len(self.temperatures))
        s_fit = np.zeros(len(self.temperatures))
        for (i, temp) in enumerate(self.temperatures):
            if temp <= T_switch:
                cp_fit[i] = a_low[0] + a_low[1] * temp + a_low[2] * temp**2.0 + a_low[3] * temp**3.0 + a_low[4] * temp**4.0
                h_fit[i] = a_low[0] * temp + a_low[1] / 2.0 * temp**2.0 + a_low[2] / 3.0 * temp**3.0 + a_low[3] / 4.0 * temp**4.0 + a_low[4] / 5.0 * temp**5.0 + a_low[5]
                s_fit[i] = a_low[0] * np.log(temp) + a_low[1] * temp + a_low[2] / 2.0 * temp**2.0 + a_low[3] / 3.0 * temp**3.0 + a_low[4] / 4.0 * temp**4.0 + a_low[6]
            else:
                cp_fit[i] = a_high[0] + a_high[1] * temp + a_high[2] * temp**2.0 + a_high[3] * temp**3.0 + a_high[4] * temp**4.0
                h_fit[i] = a_high[0] * temp + a_high[1] / 2.0 * temp**2.0 + a_high[2] / 3.0 * temp**3.0 + a_high[3] / 4.0 * temp**4.0 + a_high[4] / 5.0 * temp**5.0 + a_high[5]
                s_fit[i] = a_high[0] * np.log(temp) + a_high[1] * temp + a_high[2] / 2.0 * temp**2.0 + a_high[3] / 3.0 * temp**3.0 + a_high[4] / 4.0 * temp**4.0 + a_high[6]

        cp_fit *= R
        h_fit *= R
        s_fit *= R

        return cp_fit, h_fit, s_fit
