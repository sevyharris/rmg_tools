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
        frequencies = [49.5, 68.6, 73.6, 102.0, 437.6, 452.9, 1596.3, 3675.6, 3787.0]
        cls.adsorbate_thermo = adsorbate_thermo.AdsorbateThermoCalc(molecular_weight, frequencies)

    def test_get_translation_thermo(self):
        """
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




if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
