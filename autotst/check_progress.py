# script to check on the n-heptane model AutoTST progress
import os
import re
import glob
import rmgpy.species
import rmgpy.chemkin

# load the model
chemkin_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/chem_annotated.inp"
dictionary_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/species_dictionary.txt"
transport_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/tran.dat"
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(
    chemkin_path,
    dictionary_path=dictionary_path,
    transport_path=transport_path
)

base_dir = "/work/westgroup/harris.se/autoscience/dft/"
species_dirs = glob.glob(os.path.join(base_dir, 'thermo', 'species*'))
rmg_thermos = glob.glob(os.path.join(base_dir, 'thermo', 'species*', 'conformers', 'arkane', 'RMG_libraries', 'thermo.py'))

reaction_dirs = glob.glob(os.path.join(base_dir, 'kinetics', 'reaction_*'))
rmg_kinetics = glob.glob(os.path.join(base_dir, 'kinetics', 'reaction_*', 'conformers', 'arkane', 'RMG_libraries', 'library.py'))

completed_thermo_indices = []
completed_kinetic_indices = []

print("Thermo completed for the following species:")
for thermo_path in rmg_thermos:
    matches = re.search('species_(.*)/c', thermo_path)
    print(int(matches[1]),"\t", species_list[int(matches[1])])
    completed_thermo_indices.append(int(matches[1]))

print("Thermo incomplete for the following species:")
for species_dir in species_dirs:
    matches = re.search('species_(.*)', species_dir)
    idx = int(matches[1])
    if idx not in completed_thermo_indices:
        print(idx,"\t", species_list[idx])


print("Kinetics completed for the following reactions:")
for kinetics_path in rmg_kinetics:
    matches = re.search('reaction_(.*)/t', kinetics_path)
    idx = int(matches[1])
    print(int(matches[1]),"\t", reaction_list[idx])
    completed_kinetic_indices.append(int(matches[1]))

print("Kinetics incomplete for the following reactions:")
for reaction_dir in reaction_dirs:
    matches = re.search('reaction_(.*)', reaction_dir)
    idx = int(matches[1])
    if idx not in completed_kinetic_indices:
        print(idx,"\t", reaction_list[idx])







