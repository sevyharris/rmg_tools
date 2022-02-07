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


# print(f'R={R}')
# print(f'kB={kB}')
# print(f'h={h}')
# print(f'c={c}')
# print(f'amu={amu}')
# print(f'Avogadro={Avogadro}')


class AdsorbateThermoCalc:
    """
    Class for molecule stat thermo

    site area - assumes Pt(111) fcc facet
    """
    def __init__(
        self,
        molecular_weight,           # g/mol
        frequencies,                # TODO make a load function, pretty sure this needs to be sorted
        site_area=None,       # surface area per binding site m^2, default is Pt(111) facet
        site_occupation_number=1,   # number of sites occupied by adsorbate
        cutoff_frequency = 100.0,   # cm^-1
        twoD_gas = False,
        temperatures = None,
    ):
        self.molecular_weight = molecular_weight
        self.frequencies = frequencies
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
            dT = 10.0 #temperature increment
            temperatures = np.append(temperatures, np.arange(T_low, T_high + dT, dT))
        self.temperatures = temperatures

    
    def load_molecule_info(self):
        pass

    def get_translation_thermo(self):
        # Just using the area is not a function of temperature equations
        Q_trans = ((2.0 * pi * self.molecular_weight * amu * kB * self.temperatures) / 
                   np.float_power(h, 2.0)) * self.site_area * self.site_occupation_number
        S_trans = R * (2.0 + np.log(Q_trans))
        Cp_trans = R * np.ones(len(self.temperatures))  # NOTE: Cp = Cv 
        dH_trans = R * 1.0 * self.temperatures
        return Q_trans, S_trans, Cp_trans, dH_trans

    def get_vibrational_thermo(self):
        units = 1.0
        # TODO check units
        
        frequencies = self.frequencies  # TODO check sorted
        if self.twoD_gas:  # skip the first two if we do 2D gas
            frequencies = frequencies[2:]

        # TODO point to a LaTeX notebook explaining the equations AND the matrixization
        # x is a matrix (Temperatures * Frequencies in size)
        # x = nu * units / temp #cm^-1 * K cm / K = dimensionless
        x = np.matmul(np.matrix(units / self.temperatures).transpose(), np.matrix(frequencies))
        Q_vib = np.prod(1.0 / (1.0 - np.exp(-x)), 1)
        S_vib = np.sum(-np.log(1.0 - np.exp(-x)) + np.multiply(x, np.exp(-x)) / (1.0 - np.exp(-x)), 1) * R
        dH_vib = np.multiply(np.sum(np.multiply(x, np.exp(-x)) / (1.0 - np.exp(-x)), 1), self.temperatures) * R 
        Cv_vib = np.sum(np.multiply(np.float_power(x, 2.0), np.exp(-x)) / np.float_power((1.0 - np.exp(-x)), 2.0), 1) * R
        return Q_vib, S_vib, dH_vib, Cv_vib


def get_translation_thermo(temperatures):
    # unpack the constants (not essential, but makes it easier to read)
    # R = molecule.R
    # kB = molecule.kB
    # h = molecule.h
    # amu = molecule.amu
    # P_ref = molecule.P_ref
    # m = molecule.adsorbate_mass
    # pi = np.pi
    # area = molecule.unit_cell_area
    # sites = molecule.site_occupation_number

    Q_trans  = np.ones(len(temperatures)) 

    # Q_trans[i] = (2*pi*m*amu*kB*T/h**2) * site_area * molecule_dentate
    # S_trans[i] = R * (2.0 + np.log( Q_trans[i] ))
    # Cp_trans[i] = R * 1.0 #NOTE: Cp = Cv 
    # dH_trans[i] = R * 1.0 * T




    # # add the results to the thermo object
    # molecule.Q_trans = Q_trans
    # molecule.S_trans = S_trans
    # molecule.dH_trans = dH_trans
    # molecule.Cp_trans = Cp_trans 

