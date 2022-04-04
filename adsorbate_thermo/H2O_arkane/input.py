#!/usr/bin/env python

modelChemistry = "M06-2X/cc-pVTZ"
useHinderedRotors = True
useBondCorrections = False

frequencyScaleFactor = 0.982
species('H2O', 'conformer_0000.py', structure=SMILES('O'))

thermo('H2O', 'NASA')
