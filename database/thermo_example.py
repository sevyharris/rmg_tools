# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import rmgpy.species
import rmgpy.kinetics
import rmgpy.thermo
import rmgpy.data.rmg

import numpy as np

import matplotlib.pyplot as plt
# %matplotlib inline

# +

# Define Propane as a Species object
# -

CCC = rmgpy.species.Species(smiles='CCC')
display(CCC)

CCC.thermo

H2 = rmgpy.species.Species(smiles='[H][H]')

# pickle
# pass in arguments with sys
# pointers/memory
# chemkin files


# +
# Read in the RMG-database

# +
database = rmgpy.data.rmg.RMGDatabase()

thermo_libraries = [
    'Klippenstein_Glarborg2016',
    'BurkeH2O2',
    'thermo_DFT_CCSDTF12_BAC', 
    'DFT_QCI_thermo',
    'primaryThermoLibrary',
    'CurranPentane'
]

database.load(
    path = rmgpy.settings['database.directory'],
    thermo_libraries = thermo_libraries,
    transport_libraries = [],
    reaction_libraries = [],
    seed_mechanisms = [],
    kinetics_families = ['H_Abstraction'],
    kinetics_depositories = ['training'],
    depository = False,
)
# !!!!!!!!!!!!!!!!!!!!!!!! Keep this for kinetics families
for family in database.kinetics.families:
    if not database.kinetics.families[family].auto_generated:
        database.kinetics.families[family].add_rules_from_training(thermo_database=database.thermo)
        database.kinetics.families[family].fill_rules_by_averaging_up(verbose=True)


# -

print(database)

print(database.thermo)

print(database.kinetics)

print(database.thermo.libraries)

database.thermo.libraries['CurranPentane']

for i, key in enumerate(database.thermo.libraries['CurranPentane'].entries):
    entry = database.thermo.libraries['CurranPentane'].entries[key]
    print(i, entry)
    if i >= 1:
        break
    

print(entry.index)
print(entry.label)
print(entry.item)

print(entry.item.to_adjacency_list())

my_thermo_data = entry.data
print(my_thermo_data)

my_thermo_data.get_enthalpy(1000)

my_thermo_data.polynomials

# +
# Plot a thermo example
N = 101

Ts = np.linspace(300, 3000, N)
H = np.zeros_like(Ts)
S = np.zeros_like(Ts)
Cp = np.zeros_like(Ts)

for i in range(len(Ts)):
    H[i] = my_thermo_data.get_enthalpy(Ts[i])
    S[i] = my_thermo_data.get_entropy(Ts[i])
    Cp[i] = my_thermo_data.get_heat_capacity(Ts[i])

# -

print(my_thermo_data.get_enthalpy.__doc__)

plt.plot(Ts, H, marker='x')
plt.xlabel('T (K)')
plt.ylabel('Enthalpy (J/mol)')

my_thermo_data.get_enthalpy(298)

# +
# Get the thermo of CCC using CurranPentane library
# -

database.thermo.get_thermo_data_from_library(CCC, database.thermo.libraries['CurranPentane'])

# +
curran_H2_thermo, library, entry= database.thermo.get_thermo_data_from_library(H2, database.thermo.libraries['CurranPentane'])
burke_H2_thermo, library, entry= database.thermo.get_thermo_data_from_library(H2, database.thermo.libraries['BurkeH2O2'])

burke_H2_thermo = burke_H2_thermo.to_nasa(300, 3000, 1000)

print(curran_H2_thermo)
print()
print(burke_H2_thermo)

# +
# Compare the two
N = 101

Ts = np.linspace(300, 3000, N)
H_Curran = np.zeros_like(Ts)
H_Burke = np.zeros_like(Ts)

for i in range(len(Ts)):
    H_Curran[i] = curran_H2_thermo.get_enthalpy(Ts[i])
    H_Burke[i] = burke_H2_thermo.get_enthalpy(Ts[i])

plt.plot(Ts, H_Curran, label='Curran')
plt.plot(Ts, H_Burke, label='Burke')
plt.xlabel('T (K)')
plt.ylabel('Enthalpy (J/mol)')
plt.legend()
# -



library

entry.item

entry.label

entry.data





dir(my_thermo_data)


