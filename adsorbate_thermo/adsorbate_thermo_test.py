import os
import unittest
import numpy as np
import adsorbate_thermo

################################################################################


class AdsorbateThermoTest(unittest.TestCase):
    """
    Contains unit tests for the gaussian module, used for parsing Gaussian log files.
    """
    @classmethod
    def setUpClass(cls):
        """
        A method that is run before all unit tests in this class.
        """
        molecular_weight = 18.02  # H2O
        composition = {'H': 2, 'O': 1, 'C': 0, 'N': 0}
        frequencies = [49.5, 68.6, 73.6, 102.0, 437.6, 452.9, 1596.3, 3675.6, 3787.0]
        heat_of_formation_0K = -295.05
        cls.adsorbate_thermo = adsorbate_thermo.AdsorbateThermoCalc(molecular_weight, frequencies, composition, heat_of_formation_0K, twoD_gas=True)

    def test_get_translation_thermo(self):
        """
        Check values of Q_trans, S_trans, Cp_trans, dH_trans
        """
        Q_trans, S_trans, Cp_trans, dH_trans = self.adsorbate_thermo.get_translation_thermo()
        self.assertAlmostEqual(Q_trans[0], 121.638576209145)
        self.assertAlmostEqual(Q_trans[45], 301.9035599354932)
        self.assertAlmostEqual(Q_trans[-2], 811.8757895562588)

        self.assertAlmostEqual(S_trans[0], 56.54717436538365)
        self.assertAlmostEqual(S_trans[45], 64.10547399343213)
        self.assertAlmostEqual(S_trans[-2], 72.33048004244537)

        self.assertAlmostEqual(Cp_trans[0], 8.314472)
        self.assertAlmostEqual(Cp_trans[45], 8.314472)
        self.assertAlmostEqual(Cp_trans[-2], 8.314472)

        self.assertAlmostEqual(dH_trans[0], 2478.9598268)
        self.assertAlmostEqual(dH_trans[45], 6152.70928)
        self.assertAlmostEqual(dH_trans[-2], 16545.79928)

    def test_get_vibration_thermo(self):
        """
        Check values of Q_vib, S_vib, dH_vib, Cv_vib
        """
        Q_vib, S_vib, dH_vib, Cv_vib = self.adsorbate_thermo.get_vibrational_thermo()
        self.assertAlmostEqual(Q_vib[0], 11.034963254025913)
        self.assertAlmostEqual(Q_vib[45], 130.32233489549847)
        self.assertAlmostEqual(Q_vib[-2], 6018.728683338621)

        self.assertAlmostEqual(S_vib[0], 38.07225504419116)
        self.assertAlmostEqual(S_vib[45], 67.49815326008806)
        self.assertAlmostEqual(S_vib[-2], 110.22104390942083)

        self.assertAlmostEqual(dH_vib[0], 5399.089968555581)
        self.assertAlmostEqual(dH_vib[45], 19984.872272678833)
        self.assertAlmostEqual(dH_vib[-2], 75347.8861223393)

        self.assertAlmostEqual(Cv_vib[0], 28.074698286767678)
        self.assertAlmostEqual(Cv_vib[45], 36.77745028228141)
        self.assertAlmostEqual(Cv_vib[-2], 49.93125232130692)

    def test_get_thermo(self):
        """
        Check values of Q_vib, S_vib, dH_vib, Cv_vib
        """
        pass


if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
