#!/usr/bin/env python
# encoding: utf-8

# Define the level of theory ("model chemistry"):
modelChemistry = LevelOfTheory(method='CCSD(T)-F12', basis='cc-pVTZ-F12', software='molpro')
useHinderedRotors = False

# Define the species and TSs:
species('CH4', 'ch4.py', structure=SMILES('C'))
species('NH2', 'nh2.py', structure=SMILES('[NH2]'))
species('CH3', 'ch3.py', structure=SMILES('[CH3]'))
species('NH3', 'nh3.py', structure=SMILES('N'))
transitionState('TS1', 'ts1.py')


# Define the reactions:
reaction(
    label = 'CH4 + NH2 <=> CH3 + NH3',
    reactants = ['CH4', 'NH2'],
    products = ['CH3', 'NH3'],
    transitionState = 'TS1',
#    tunneling='Eckart',
)


# Request rate coefficient computations:
kinetics(label='CH4 + NH2 <=> CH3 + NH3',
         Tmin=(300,'K'), Tmax=(2500,'K'), Tcount = 25)

