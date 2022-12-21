import os
import copy
import itertools

from rmgpy.exceptions import ActionError
import rmgpy.reaction
import rmgpy.chemkin
import rmgpy.data.kinetics


# Generated using get_reactant_and_product_numbers.py
reactant_structures = {'1,2_XY_interchange': [1, 1], '1,2_shiftS': [1, 1], '1,4_Linear_birad_scission': [1, 2], 'Cyclic_Thioether_Formation': [1, 2], 'Intra_RH_Add_Endocyclic': [1, 1], 'Intra_RH_Add_Exocyclic': [1, 1], 'Intra_R_Add_ExoTetCyclic': [1, 2], 'Intra_Retro_Diels_alder_bicyclic': [1, 1], 'Korcek_step1': [1, 1], 'Korcek_step1_cat': [2, 2], 'Korcek_step2': [1, 2], 'Surface_DoubleBond_to_Bidentate': [2, 1], 'intra_substitutionCS_cyclization': [1, 1], 'lone_electron_pair_bond': [1, 1], '1+2_Cycloaddition': [2, 1], '1,2-Birad_to_alkene': [1, 1], '1,2_Insertion_CO': [2, 1], '1,2_Insertion_carbene': [2, 1], '1,2_NH3_elimination': [1, 2], '1,2_shiftC': [1, 1], '1,3_Insertion_CO2': [2, 1], '1,3_Insertion_ROR': [2, 1], '1,3_Insertion_RSR': [2, 1], '1,3_NH3_elimination': [1, 2], '1,3_sigmatropic_rearrangement': [1, 1], '1,4_Cyclic_birad_scission': [1, 1], '2+2_cycloaddition': [2, 1], '6_membered_central_C-C_shift': [1, 1], 'Baeyer-Villiger_step1_cat': [3, 2], 'Baeyer-Villiger_step2': [1, 2], 'Baeyer-Villiger_step2_cat': [2, 3], 'Bimolec_Hydroperoxide_Decomposition': [2, 3], 'Birad_R_Recombination': [2, 1], 'Birad_recombination': [1, 1], 'Br_Abstraction': [2, 2], 'CO_Disproportionation': [2, 2], 'Cl_Abstraction': [2, 2], 'Concerted_Intra_Diels_alder_monocyclic_1,2_shiftH': [1, 1], 'Cyclic_Ether_Formation': [1, 2], 'Cyclopentadiene_scission': [1, 1], 'Diels_alder_addition': [2, 1], 'Diels_alder_addition_Aromatic': [2, 1], 'Disproportionation': [2, 2], 'Disproportionation-Y': [2, 2], 'F_Abstraction': [2, 2], 'H2_Loss': [1, 2], 'HO2_Elimination_from_PeroxyRadical': [1, 2], 'H_Abstraction': [2, 2], 'Intra_2+2_cycloaddition_Cd': [1, 1], 'Intra_5_membered_conjugated_C=C_C=C_addition': [1, 1], 'Intra_Diels_alder_monocyclic': [1, 1], 'Intra_Disproportionation': [1, 1], 'Intra_R_Add_Endocyclic': [1, 1], 'Intra_R_Add_Exo_scission': [1, 1], 'Intra_R_Add_Exocyclic': [1, 1], 'Intra_ene_reaction': [1, 1], 'Ketoenol': [1, 1], 'Peroxyl_Disproportionation': [2, 3], 'Peroxyl_Termination': [2, 3], 'R_Addition_COm': [2, 1], 'R_Addition_CSm': [2, 1], 'R_Addition_MultipleBond': [2, 1], 'R_Recombination': [2, 1], 'Retroene': [1, 2], 'Singlet_Carbene_Intra_Disproportionation': [1, 1], 'Singlet_Val6_to_triplet': [1, 1], 'SubstitutionS': [2, 2], 'Substitution_O': [2, 2], 'Surface_Abstraction': [2, 2], 'Surface_Abstraction_Beta': [2, 2], 'Surface_Abstraction_Beta_double_vdW': [2, 2], 'Surface_Abstraction_Single_vdW': [2, 2], 'Surface_Abstraction_vdW': [2, 2], 'Surface_Addition_Single_vdW': [2, 2], 'Surface_Adsorption_Abstraction_vdW': [2, 2], 'Surface_Adsorption_Bidentate': [3, 1], 'Surface_Adsorption_Dissociative': [3, 2], 'Surface_Adsorption_Dissociative_Double': [3, 2], 'Surface_Adsorption_Double': [2, 1], 'Surface_Adsorption_Single': [2, 1], 'Surface_Adsorption_vdW': [2, 1], 'Surface_Bidentate_Dissociation': [2, 1], 'Surface_Dissociation': [2, 2], 'Surface_Dissociation_Beta': [2, 2], 'Surface_Dissociation_Double_vdW': [2, 2], 'Surface_Dissociation_to_Bidentate': [3, 2], 'Surface_Dissociation_vdW': [2, 2], 'Surface_Dual_Adsorption_vdW': [2, 2], 'Surface_EleyRideal_Addition_Multiple_Bond': [2, 1], 'Surface_Migration': [1, 1], 'Surface_vdW_to_Bidentate': [2, 1], 'XY_Addition_MultipleBond': [2, 1], 'XY_elimination_hydroxyl': [1, 3], 'halocarbene_recombination': [2, 1], 'halocarbene_recombination_double': [2, 1], 'intra_H_migration': [1, 1], 'intra_NO2_ONO_conversion': [1, 1], 'intra_OH_migration': [1, 1], 'intra_halogen_migration': [1, 1], 'intra_substitutionCS_isomerization': [1, 1], 'intra_substitutionS_cyclization': [1, 2], 'intra_substitutionS_isomerization': [1, 1]}


def relabel2_2(input_r, family):
    input_reaction = copy.deepcopy(input_r)

    # copied from AutoTST.autotst.reaction.py
    def get_rmg_mol(smile):
        smiles_conversions = {
            "[CH]": "[CH...]",
            "CARBONMONOXIDE": "[C-]#[O+]"
        }

        if smile.upper() in list(smiles_conversions.keys()):
            smile = smiles_conversions[smile.upper()]
        return rmgpy.molecule.Molecule(smiles=smile).generate_resonance_structures()

    rmg_reactants = [get_rmg_mol(sp.smiles) for sp in input_reaction.reactants]
    rmg_products = [get_rmg_mol(sp.smiles) for sp in input_reaction.products]

    combos_to_try = list(itertools.product(
        list(itertools.product(*rmg_reactants)),
        list(itertools.product(*rmg_products))
    ))

    for rmg_reactants, rmg_products in combos_to_try:

        test_reaction = rmgpy.reaction.Reaction(
            reactants=list(rmg_reactants),
            products=list(rmg_products)
        )

        try:
            labeled_r, labeled_p = ref_db.kinetics.families[family].get_labeled_reactants_and_products(
                test_reaction.reactants,
                test_reaction.products
            )

            if labeled_r is None or labeled_p is None:
                continue

            if input_reaction.reactants[0].molecule[0].is_isomorphic(labeled_r[0]):
                input_reaction.reactants[0].molecule[0] = labeled_r[0]
                input_reaction.reactants[1].molecule[0] = labeled_r[1]
            else:
                input_reaction.reactants[0].molecule[0] = labeled_r[1]
                input_reaction.reactants[1].molecule[0] = labeled_r[0]

            if input_reaction.products[0].molecule[0].is_isomorphic(labeled_p[0]):
                input_reaction.products[0].molecule[0] = labeled_p[0]
                input_reaction.products[1].molecule[0] = labeled_p[1]
            else:
                input_reaction.products[0].molecule[0] = labeled_p[1]
                input_reaction.products[1].molecule[0] = labeled_p[0]

            return input_reaction

        except ActionError:
            pass
    return False


def get_family(reaction):
    n_reactants = len(reaction.reactants)
    n_products = len(reaction.products)
    print('possible families:')
    for family in reactant_structures.keys():
        if reactant_structures[family][0] == n_reactants and reactant_structures[family][1] == n_products:
            print(family)

    return

    include_list = [  # these are bimoleculars
        'Disproportionation',
        'Disproportionation-Y',
        'F_Abstraction',
        'H_Abstraction', 
    ]
    for family in ref_db.kinetics.families:
#         if family in skip_list:
#             continue
        if family not in include_list:
            continue
            
            
#         print(f'Trying family: {family}')  # for debugging
        my_family = ref_db.kinetics.families[family]

        test_rxn = relabel(reaction, family)
        if test_rxn is False:
            continue

#         # print atom labels to check labeling
#         for sp in test_rxn.reactants + test_rxn.products:
#             print(sp.molecule[0].get_all_labeled_atoms())

        try:
            template = my_family.get_reaction_template(test_rxn)
            return family
        except rmgpy.exceptions.UndeterminableKineticsError:
            continue


if __name__ == '__main__':
    # Run example

    # load minimal ethane example
    chemkin_file = '/home/moon/rmg/my_examples/ethane/chemkin/chem_annotated.inp'
    dict_file = '/home/moon/rmg/my_examples/ethane/chemkin/species_dictionary.txt'
    species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dict_file, use_chemkin_names=True)

    ###############################################################
    # load the database
    # load the thermo database,
    # this is minimal for better speed
    thermo_libs = [
        'primaryThermoLibrary',
    ]

    thermo_library_path = os.path.join(rmgpy.settings['database.directory'], 'thermo')
    thermo_database = rmgpy.data.thermo.ThermoDatabase()
    thermo_database.load(
        thermo_library_path,
        libraries=thermo_libs
    )

    # load the families
    ref_library_path = os.path.join(rmgpy.settings['database.directory'], 'kinetics')
    kinetics_database = rmgpy.data.kinetics.KineticsDatabase()
    kinetics_database.load(
        ref_library_path,
        libraries=[],
        families='all'
    )

    # load the entire database
    ref_db = rmgpy.data.rmg.RMGDatabase()
    ref_db.kinetics = kinetics_database
    ref_db.thermo = thermo_database

    for family in ref_db.kinetics.families:
        my_family = ref_db.kinetics.families[family]

        if my_family.product_num is None or my_family.reactant_num is None:
            print(my_family)

    #################################################################
    # construct a reaction by hand (or load from the NIST model)
    reactant1 = rmgpy.species.Species(smiles='[CH3]')
    reactant2 = rmgpy.species.Species(smiles='CC')
    product1 = rmgpy.species.Species(smiles='C')
    product2 = rmgpy.species.Species(smiles='C[CH2]')

    reactant1.thermo = ref_db.thermo.get_thermo_data(reactant1)
    reactant2.thermo = ref_db.thermo.get_thermo_data(reactant2)
    product1.thermo = ref_db.thermo.get_thermo_data(product1)
    product2.thermo = ref_db.thermo.get_thermo_data(product2)

    my_reaction = rmgpy.reaction.Reaction()
    my_reaction.reactants = [reactant1, reactant2]
    my_reaction.products = [product1, product2]

    # # reverse it
    # my_reaction.products = [reactant1, reactant2]
    # my_reaction.reactants = [product1, product2]

    #################################################################
    family = get_family(reaction_list[1])
    print(family)
