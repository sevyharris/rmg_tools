# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import os
import subprocess
import pickle
import cantera as ct
import rmgpy.chemkin
import numpy as np


# +
def names_match(rmg_rxn, ct_rxn):
    rmg_r = set([str(x) for x in rmg_rxn.reactants])
    rmg_p = set([str(x) for x in rmg_rxn.products])

    ct_r = set(ct_rxn.reactants.keys())
    ct_p = set(ct_rxn.products.keys())

    return rmg_r == ct_r and rmg_p == ct_p

def ct_arrhenius_match(ct_rate1, ct_rate2):
    # assert isinstance(ct_rate1, ct.reaction.ArrheniusRate)
    # assert isinstance(ct_rate2, ct.reaction.ArrheniusRate)
    if not np.isclose(ct_rate1.pre_exponential_factor, ct_rate2.pre_exponential_factor):
        return False
    if not np.isclose(ct_rate1.temperature_exponent, ct_rate2.temperature_exponent):
        return False
    if not np.isclose(ct_rate1.activation_energy, ct_rate2.activation_energy):
        return False
    return True



def kinetics_match(rmg_rxn, ct_rxn):
    if isinstance(rmg_rxn.kinetics, rmgpy.kinetics.arrhenius.MultiArrhenius) and isinstance(ct_rxn.rate, ct.reaction.ArrheniusRate):
        for i in range(len(rmg_rxn.kinetics.arrhenius)):
            if ct_arrhenius_match(rmg_rxn.kinetics.arrhenius[i].to_cantera_kinetics(), ct_rxn.rate):
                return True
        return False
    if isinstance(rmg_rxn.kinetics, rmgpy.kinetics.arrhenius.MultiPDepArrhenius) and isinstance(ct_rxn.rate, ct.reaction.PlogRate):
        # if any of the PLOGs match, it's a match
        for i in range(len(rmg_rxn.kinetics.arrhenius)):
            dummy_rmg_rxn = rmgpy.reaction.Reaction()
            dummy_rmg_rxn.kinetics = rmg_rxn.kinetics.arrhenius[i]
            if kinetics_match(dummy_rmg_rxn, ct_rxn):
                return True
        return False
    if isinstance(rmg_rxn.kinetics, rmgpy.kinetics.arrhenius.Arrhenius) and isinstance(ct_rxn.rate, ct.reaction.ArrheniusRate):
        return ct_arrhenius_match(rmg_rxn.kinetics.to_cantera_kinetics(), ct_rxn.rate)
    if isinstance(rmg_rxn.kinetics, rmgpy.kinetics.falloff.ThirdBody) and isinstance(ct_rxn.rate, ct.reaction.ArrheniusRate):
        return ct_arrhenius_match(rmg_rxn.kinetics.arrheniusLow.to_cantera_kinetics(), ct_rxn.rate)
    if isinstance(rmg_rxn.kinetics, rmgpy.kinetics.falloff.Troe) and isinstance(ct_rxn.rate, ct.reaction.TroeRate):
        if not ct_arrhenius_match(rmg_rxn.kinetics.arrheniusLow.to_cantera_kinetics(), ct_rxn.rate.low_rate):
            return False
        if not np.isclose(rmg_rxn.kinetics.alpha, ct_rxn.rate.falloff_coeffs[0]):
            return False
        return True
    if isinstance(rmg_rxn.kinetics, rmgpy.kinetics.chebyshev.Chebyshev) and isinstance(ct_rxn.rate, ct.reaction.ChebyshevRate):
        return np.sum(np.isclose(rmg_rxn.kinetics.coeffs.value_si, ct_rxn.rate.data, rtol=1e-2, atol=1e-3)) > 8 # TODO figure out why the Cantera/RMG discrepancy for RMG 38, Cantera 47
    if isinstance(rmg_rxn.kinetics, rmgpy.kinetics.arrhenius.PDepArrhenius) and isinstance(ct_rxn.rate, ct.reaction.PlogRate):
        # have to collect all arrheniuses, including MultiArrheniuses
        rmg_arrheniuses = []
        rmg_pressures = []
        for i in range(len(rmg_rxn.kinetics.arrhenius)):
            if isinstance(rmg_rxn.kinetics.arrhenius[i], rmgpy.kinetics.arrhenius.MultiArrhenius):
                for j in range(len(rmg_rxn.kinetics.arrhenius[i].arrhenius)):
                    rmg_pressures.append(rmg_rxn.kinetics.pressures.value_si[i])
                    rmg_arrheniuses.append(rmg_rxn.kinetics.arrhenius[i].arrhenius[j])

            else:
                rmg_arrheniuses.append(rmg_rxn.kinetics.arrhenius[i])
                rmg_pressures.append(rmg_rxn.kinetics.pressures.value_si[i])
        if len(rmg_arrheniuses) != len(ct_rxn.rate.rates):
            return False
        for i in range(len(rmg_arrheniuses)):
            if not np.isclose(rmg_pressures[i], ct_rxn.rate.rates[i][0]):
                return False
            if not ct_arrhenius_match(rmg_arrheniuses[i].to_cantera_kinetics(), ct_rxn.rate.rates[i][1]):
                return False
        return True
    if isinstance(rmg_rxn.kinetics, rmgpy.kinetics.falloff.Lindemann) and isinstance(ct_rxn.rate, ct.reaction.LindemannRate):
        if not ct_arrhenius_match(rmg_rxn.kinetics.arrheniusLow.to_cantera_kinetics(), ct_rxn.rate.low_rate):
            return False
        return ct_arrhenius_match(rmg_rxn.kinetics.arrheniusHigh.to_cantera_kinetics(), ct_rxn.rate.high_rate)
    return False


# -

kinetics_match(reaction_list[125], gas.reactions()[134])

for i in range(len(reaction_list)):
    if isinstance(reaction_list[i].kinetics, rmgpy.kinetics.arrhenius.MultiPDepArrhenius):
        print(i)
        break

reaction_list[125].kinetics.arrhenius[0]

rmg2ct_possibilities[125]

gas.reactions()[134].rate.rates


# +
def rmg2str(rmg_reaction):
    prod_str = '+'.join(sorted([x.to_chemkin() for x in rmg_reaction.products]))
    react_str = '+'.join(sorted([x.to_chemkin() for x in rmg_reaction.reactants]))

    rxn_str = '<=>'.join(sorted([react_str, prod_str]))
    return rxn_str

def rmg2str_cancel3rd(rmg_reaction):
    # cancel species that appear on both sides of the reactions
    products = sorted([x.to_chemkin() for x in rmg_reaction.products])
    reactants = sorted([x.to_chemkin() for x in rmg_reaction.reactants])

    items_to_remove = []
    for r in reactants:
        if r in products:
            items_to_remove.append(r)
            break
    if items_to_remove:
        products.remove(items_to_remove[0])
        reactants.remove(items_to_remove[0])
    prod_str = '+'.join(products)
    react_str = '+'.join(reactants)

    rxn_str = '<=>'.join(sorted([react_str, prod_str]))
    return rxn_str

def ct2str(ct_reaction):
    r_list = []
    for reactant in ct_reaction.reactants.keys():
        for j in range(int(ct_reaction.reactants[reactant])):
            r_list.append(reactant)
    p_list = []
    for product in ct_reaction.products.keys():
        for j in range(int(ct_reaction.products[product])):
            p_list.append(product)
    prod_str = '+'.join(sorted([p for p in p_list]))
    react_str = '+'.join(sorted([r for r in r_list]))

    rxn_str = '<=>'.join(sorted([react_str, prod_str]))
    return rxn_str



# +
chemkin_file = '/home/moon/uq_reformulation_paper/01_RMG_propane_model/chem_annotated.inp'
species_dictionary = '/home/moon/uq_reformulation_paper/01_RMG_propane_model/species_dictionary.txt'
cantera_file = '/home/moon/uq_reformulation_paper/01_RMG_propane_model/chem_annotated.yaml'
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, species_dictionary)
gas = ct.Solution(cantera_file)


# -

rmg_rxn_strs = [rmg2str(rxn) for rxn in reaction_list]
rmg_rxn_strs_cancelled = [rmg2str_cancel3rd(rxn) for rxn in reaction_list]
ct_rxn_strs = [ct2str(gas.reactions()[i]) for i in range(gas.n_reactions)]

# +
rmg2ct_possibilities = []

for i in range(len(reaction_list)):
    index_possible = []
    for j in range(gas.n_reactions):
        if rmg_rxn_strs[i] == ct_rxn_strs[j]:
            index_possible.append(j)
    rmg2ct_possibilities.append(index_possible)

    if not index_possible:
        print(i, reaction_list[i])
        print()
# -

ct2rmg_mapping = {}
for i in range(len(rmg2ct_possibilities)):
    if len(rmg2ct_possibilities[i]) == 1:
        assert kinetics_match(reaction_list[i], gas.reactions()[rmg2ct_possibilities[i][0]]), i
        # TODO check that the kinetics match too
        ct2rmg_mapping[rmg2ct_possibilities[i][0]] = i

remaining_ct = set(np.arange(gas.n_reactions)) - set(ct2rmg_mapping.keys())
remaining_rmg = set(np.arange(len(reaction_list))) - set(ct2rmg_mapping.values())

# +
# try to match the rest based on kinetics
# -

# # Solve the ones that have zero possibilities

# +
for i in remaining_rmg:
    if rmg2ct_possibilities[i]:
        continue
    new_possibilities = []
    for j in range(gas.n_reactions):
        if rmg_rxn_strs_cancelled[i] == ct_rxn_strs[j] and kinetics_match(reaction_list[i], gas.reactions()[j]):
            # print(i, j)
            new_possibilities.append(j)
    if len(new_possibilities) == 1:
        ct2rmg_mapping[new_possibilities[0]] = i

remaining_ct = set(np.arange(gas.n_reactions)) - set(ct2rmg_mapping.keys())
remaining_rmg = set(np.arange(len(reaction_list))) - set(ct2rmg_mapping.values())
# -

# # Handle multiarrhenius

# +
for i in remaining_rmg:
    if not isinstance(reaction_list[i].kinetics, rmgpy.kinetics.arrhenius.MultiArrhenius):
        continue
    possibilities = []
    for j in remaining_ct:
        if kinetics_match(reaction_list[i], gas.reactions()[j]):
            possibilities.append(j)
    if len(possibilities) == len(reaction_list[i].kinetics.arrhenius):
        for j in possibilities:
            ct2rmg_mapping[j] = i

remaining_ct = set(np.arange(gas.n_reactions)) - set(ct2rmg_mapping.keys())
remaining_rmg = set(np.arange(len(reaction_list))) - set(ct2rmg_mapping.values())
# -



# # Handle multiarrhenius PDEP

# +
for i in remaining_rmg:
    if not isinstance(reaction_list[i].kinetics, rmgpy.kinetics.arrhenius.MultiPDepArrhenius):
        continue
    possibilities = []
    for j in remaining_ct:
        if kinetics_match(reaction_list[i], gas.reactions()[j]):
            possibilities.append(j)
    if len(possibilities) == len(reaction_list[i].kinetics.arrhenius):
        for j in possibilities:
            ct2rmg_mapping[j] = i

remaining_ct = set(np.arange(gas.n_reactions)) - set(ct2rmg_mapping.keys())
remaining_rmg = set(np.arange(len(reaction_list))) - set(ct2rmg_mapping.values())
# -

# # look for plain kinetics matches among the survivors

for i in remaining_rmg:
    possibilities = []
    for j in remaining_ct:
        if kinetics_match(reaction_list[i], gas.reactions()[j]):
            possibilities.append(j)
    if len(possibilities) == 1:
        ct2rmg_mapping[possibilities[0]] = i
remaining_ct = set(np.arange(gas.n_reactions)) - set(ct2rmg_mapping.keys())
remaining_rmg = set(np.arange(len(reaction_list))) - set(ct2rmg_mapping.values())

print(len(remaining_ct))
print(len(remaining_rmg))

# +
# do it again? see if we've eliminated anything
# -

for i in remaining_rmg:
    possibilities = []
    for j in remaining_ct:
        if kinetics_match(reaction_list[i], gas.reactions()[j]) and (rmg_rxn_strs_cancelled[i] == ct_rxn_strs[j]):
            possibilities.append(j)
    if len(possibilities) == 1:
        ct2rmg_mapping[possibilities[0]] = i
remaining_ct = set(np.arange(gas.n_reactions)) - set(ct2rmg_mapping.keys())
remaining_rmg = set(np.arange(len(reaction_list))) - set(ct2rmg_mapping.values())

print(len(remaining_ct))
print(len(remaining_rmg))

# # Take care of Ar/He distinction

for i in remaining_rmg:
    possibilities = []
    for j in remaining_ct:
        if kinetics_match(reaction_list[i], gas.reactions()[j]) and (rmg_rxn_strs_cancelled[i] == ct_rxn_strs[j]):
            if 'He' in rmg_rxn_strs[i] and 'He' in gas.reaction_equations()[j]:
                possibilities.append(j)
            if 'Ar' in rmg_rxn_strs[i] and 'Ar' in gas.reaction_equations()[j]:
                possibilities.append(j)

    if len(possibilities) == 1:
        ct2rmg_mapping[possibilities[0]] = i
remaining_ct = set(np.arange(gas.n_reactions)) - set(ct2rmg_mapping.keys())
remaining_rmg = set(np.arange(len(reaction_list))) - set(ct2rmg_mapping.values())

print(len(remaining_ct))
print(len(remaining_rmg))

# +
# run it again
# -

for i in remaining_rmg:
    possibilities = []
    for j in remaining_ct:
        if kinetics_match(reaction_list[i], gas.reactions()[j]) and (rmg_rxn_strs_cancelled[i] == ct_rxn_strs[j]):
            possibilities.append(j)
    if len(possibilities) == 1:
        ct2rmg_mapping[possibilities[0]] = i
remaining_ct = set(np.arange(gas.n_reactions)) - set(ct2rmg_mapping.keys())
remaining_rmg = set(np.arange(len(reaction_list))) - set(ct2rmg_mapping.values())

print(len(remaining_ct))
print(len(remaining_rmg))

for i in remaining_rmg:
    print(i, type(reaction_list[i].kinetics))

reaction_list[1030].kinetics

reaction_list[870].kinetics

reaction_list[871].kinetics

#
# # At this point, just pick the thing that comes first because the order's more likely to line up

# +
for i in remaining_rmg:
    # print(i, rmg2ct_possibilities[i])

    # make the assignment
    ct2rmg_mapping[rmg2ct_possibilities[i][0]] = i
    rmg2ct_possibilities[i] = [rmg2ct_possibilities[i][0]]

    # update rmg2ct_possibilities
    for j in remaining_rmg:
        if j == i:
            continue
        try:
            rmg2ct_possibilities[j].remove(rmg2ct_possibilities[i][0])   
        except ValueError:
            pass

remaining_ct = set(np.arange(gas.n_reactions)) - set(ct2rmg_mapping.keys())
remaining_rmg = set(np.arange(len(reaction_list))) - set(ct2rmg_mapping.values())
# -

remaining_rmg

remaining_ct

for ct_index, rmg_index in ct2rmg_mapping.items():
    assert kinetics_match(reaction_list[rmg_index], gas.reactions()[ct_index])
    assert ct_rxn_strs[ct_index] == rmg_rxn_strs_cancelled[rmg_index] or ct_rxn_strs[ct_index] == rmg_rxn_strs[rmg_index]

with open(os.path.join(os.path.dirname(chemkin_file), 'ct2rmg_map.pickle'), 'wb') as f:
    pickle.dump(ct2rmg_mapping, f)

with open(os.path.join(os.path.dirname(chemkin_file), 'ct2rmg_map.pickle'), 'rb') as f:
    new_map = pickle.load(f)

assert set(new_map.keys()) == set(np.arange(gas.n_reactions))
assert set(new_map.values()) == set(np.arange(len(reaction_list)))

ct2rmg_matrix = np.zeros((gas.n_reactions, len(reaction_list)))
for ct_idx, rmg_idx in ct2rmg_rxn.items():
    ct2rmg_matrix[ct_idx, rmg_idx] = 1.0

ct2rmg_matrix.shape

reaction_list[1].kinetics.arrhenius


