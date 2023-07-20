#!/usr/bin/env python
# coding: utf-8

import os

import rmgpy.reaction
import rmgpy.data.kinetics
import rmgpy.data.thermo
import rmgpy.data.rmg
import rmgpy.molecule
import rmgpy.species



# load the thermo and kinetics databases
thermo_library_path = os.path.join(rmgpy.settings['database.directory'], 'thermo')
thermo_database = rmgpy.data.thermo.ThermoDatabase()
thermo_database.load(
    thermo_library_path,
    libraries=['primaryThermoLibrary']
)

# load the families
ref_library_path = os.path.join(rmgpy.settings['database.directory'], 'kinetics')
kinetics_database = rmgpy.data.kinetics.KineticsDatabase()
kinetics_database.load(
    ref_library_path,
    libraries=[],
    families='Disproportionation'
)

# load the entire database
ref_db = rmgpy.data.rmg.RMGDatabase()
ref_db.kinetics = kinetics_database
ref_db.thermo = thermo_database


r1 = rmgpy.species.Species(smiles='[O][O]')
r2 = rmgpy.species.Species(smiles='C=C[O]')
p1 = rmgpy.species.Species(smiles='[O]O')
p2 = rmgpy.species.Species(smiles='C=C=O')

r1.thermo = ref_db.thermo.get_thermo_data_from_groups(r1)
r2.thermo = ref_db.thermo.get_thermo_data_from_groups(r2)
p1.thermo = ref_db.thermo.get_thermo_data_from_groups(p1)
p2.thermo = ref_db.thermo.get_thermo_data_from_groups(p2)


rxn = rmgpy.reaction.Reaction()
rxn.reactants = [r1, r2]
rxn.products = [p1, p2]


rxn.reactants[0].molecule[0].atoms[0].label = '*1'
rxn.reactants[1].molecule[0].atoms[5].label = '*4'
rxn.reactants[1].molecule[0].atoms[0].label = '*3'
rxn.reactants[1].molecule[0].atoms[2].label = '*2'


rxn.reactants[0].molecule[0].get_all_labeled_atoms()


rxn.reactants[1].molecule[0].get_all_labeled_atoms()


ref_db.kinetics.families['Disproportionation'].get_reaction_template(rxn)

