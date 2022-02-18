import glob
import os
import numpy as np

import ase.io.gaussian
from ase.calculators.calculator import PropertyNotImplementedError

import rmgpy.species
import rmgpy.chemkin

import autotst.species
from autotst.calculator.gaussian import Gaussian
import autotst.calculator.statmech


# load the model
chemkin_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/chem_annotated.inp"
dictionary_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/species_dictionary.txt"
transport_path = "/home/harris.se/rmg/rmg_tools/uncertainty/nheptane/tran.dat"
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(
    chemkin_path,
    dictionary_path=dictionary_path,
    transport_path=transport_path
)

species_indices = [7, 17, 58, 68, 69, 71, 73, 77, 79, 80, 84]
for species_index in species_indices:

    conformer_dir = f"/work/westgroup/harris.se/autoscience/dft/thermo/species_{species_index:04}/conformers/"
    conformer_files = glob.glob(os.path.join(conformer_dir, 'conformer_*.log'))

    energies = np.zeros(len(conformer_files))
    for i, conformer_file in enumerate(conformer_files):
        try:
            with open(conformer_file, 'r') as f:
                atoms = ase.io.gaussian.read_gaussian_out(f)
                energy = atoms.get_potential_energy()
                energies[i] = energy
        except IndexError:
            pass
        except PropertyNotImplementedError:
            pass

    lowest_idx = np.argmin(energies)

    # read in the lowest energy conformer
    with open(conformer_files[lowest_idx], 'r') as f:
        atoms = ase.io.gaussian.read_gaussian_out(f)

    # make a conformer object again
    new_species = autotst.species.Species()
    new_cf = autotst.species.Conformer(rmg_molecule=species_list[species_index].molecule[0])
    new_cf._ase_molecule = atoms
    new_cf.update_coords_from(mol_type="ase")
    # new_cf.view()

    reaction = None
    arkane_dir = os.path.join(conformer_dir, 'arkane')
    os.makedirs(arkane_dir, exist_ok=True)
    stat = autotst.calculator.statmech.StatMech(reaction, directory=arkane_dir)
    stat.write_conformer_file2(new_cf, arkane_dir, conformer_files[lowest_idx], include_rotors=False)

    # write input file
    input_file = os.path.join(arkane_dir, 'input.py')
    formula = new_cf.rmg_molecule.get_formula()
    lines = [
        '#!/usr/bin/env python\n\n',
        f'modelChemistry = "M06-2X/cc-pVTZ"\n',
        'useHinderedRotors = False\n',
        'useBondCorrections = False\n\n',
        
        'frequencyScaleFactor = 0.982\n',

        f"species('{formula}', '{os.path.basename(conformer_files[lowest_idx][:-4])}.py', structure=SMILES('{new_cf.rmg_molecule.smiles}'))\n\n",

        f"thermo('{formula}', 'NASA')\n",
    ]
    with open(input_file, 'w') as f:
        f.writelines(lines)
        
