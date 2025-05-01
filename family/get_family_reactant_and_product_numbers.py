import os
import numpy as np
import scipy.stats
import rmgpy.chemkin
import rmgpy.data.kinetics


if __name__ == '__main__':
    reaction_structures = {  # these have no training rxns, have to be done manually
        '1,2_XY_interchange': [1, 1],
        '1,2_shiftS': [1, 1],
        '1,4_Linear_birad_scission': [1, 2],
        'Cyclic_Thioether_Formation': [1, 2],
        'Intra_RH_Add_Endocyclic': [1, 1],
        'Intra_RH_Add_Exocyclic': [1, 1],
        'Intra_R_Add_ExoTetCyclic': [1, 2],
        'Intra_Retro_Diels_alder_bicyclic': [1, 1],
        'Korcek_step1': [1, 1],
        'Korcek_step1_cat': [2, 2],
        'Korcek_step2': [1, 2],
        'Surface_DoubleBond_to_Bidentate': [2, 1],
        'intra_substitutionCS_cyclization': [1, 1],
        'lone_electron_pair_bond': [1, 1],
    }  

    # load the families
    ref_library_path = os.path.join(rmgpy.settings['database.directory'], 'kinetics')
    kinetics_database = rmgpy.data.kinetics.KineticsDatabase()
    kinetics_database.load(
        ref_library_path,
        libraries=[],
        families='all'
    )

    for family in kinetics_database.families:
        my_family = kinetics_database.families[family]

        if my_family.product_num is None or my_family.reactant_num is None:
            # print(my_family)
            # need to load training depo
            training_depo = my_family.get_training_depository()

            if len(training_depo.entries) == 0:
                # print(f'{family}')  # for manual entry
                continue

            n_reactants = []
            n_products = []
            for entry in training_depo.entries:
                rxn = training_depo.entries[entry].item
                n_reactants.append(len(rxn.reactants))
                n_products.append(len(rxn.products))

            reaction_structures[family] = [int(scipy.stats.mode(n_reactants)[0]), int(scipy.stats.mode(n_products)[0])]
            # print('\t', len(rxn.reactants), len(rxn.products))

        else:
            # print(f'{my_family} has reactant_num = {my_family.reactant_num} and product_num = {my_family.product_num}')
            reaction_structures[family] = [my_family.reactant_num, my_family.product_num]
    print(reaction_structures)
