# Script to fix incorrectly marked duplicate reactions for a chemkin file

import os
import sys


chemkin_file = sys.argv[1]

# I know it's bad form, but
# lines, reactions, and reaction_names are global variables
with open(chemkin_file, 'r') as f:
    lines = f.readlines()

reactions = []
for i, line in enumerate(lines):
    if line[0] == '!':
        continue
    trimmed_line = line.split('!')[0]  # remove anything that comes after a comment
    if '=' in trimmed_line:
        reactions.append((i, line))
    # if line[0] != '!' and '=' in line:
    #     reactions.append((i, line))

print(f'Found {len(reactions)} reactions')
# print(reactions)
reaction_names = [rxn[1].split()[0] for rxn in reactions]


def same_reaction(name1, name2):
    # name1 = name1.split()[0]
    # name1 = name1.split()[0]
    reversible1 = False
    reversible2 = False
    # TODO manage reversibility

    half_tokens1 = name1.split("=>")
    if len(half_tokens1) == 1:
        half_tokens1 = name1.split("=")
        reversible1 = True
    if len(half_tokens1) == 1:
        print('1:\t', name1)
        print('2:\t', name2)
        raise ValueError

    half_tokens2 = name2.split("=>")
    if len(half_tokens2) == 1:
        half_tokens2 = name2.split("=")
        reversible2 = True
    if len(half_tokens2) == 1:
        raise ValueError

    # TODO case where one says A + A and the other says 2A

    reactants_tokens1 = set(half_tokens1[0].split('+'))
    reactants_tokens2 = set(half_tokens2[0].split('+'))

    products_tokens1 = set(half_tokens1[1].split('+'))
    products_tokens2 = set(half_tokens2[1].split('+'))

    for r in reactants_tokens1:
        if r[0] == '2':
            raise NotImplementedError

    duplicate = False
    if reactants_tokens1 == reactants_tokens2 and products_tokens1 == products_tokens2:
        duplicate = True
    elif reversible1 and reversible2 and \
            reactants_tokens1 == products_tokens2 and products_tokens1 == reactants_tokens2:
        duplicate = True

    return duplicate


duplicate_names = set()
duplicates = []
for i in range(0, len(reaction_names)):
    for j in range(0, i):
        # if reaction_names[i] == reaction_names[j]:
        if same_reaction(reaction_names[i], reaction_names[j]):
            duplicates.append((j, i))
            #print('Found duplicates:')
            #print('\t', reaction_names[i])
            #print('\t', reaction_names[j])
            #print()
            duplicate_names.add(reaction_names[j])
            duplicate_names.add(reaction_names[i])


def marked_duplicate(reaction_number):
    rxn1_line = reactions[reaction_number][0]
    for a in range(rxn1_line, min(rxn1_line + 20, len(lines))):
        if lines[a].strip() == 'DUPLICATE':
            return True
        elif a > rxn1_line and '=' in lines[a]:
            return False
        # elif lines[a].strip() == '':
        #     return False
        # elif lines[a][0] != '!' and '=' in lines[a][0]:
        #     return False
    return False

    # rxn1_line = reactions[reaction_number][0]
    # for a in range(rxn1_line, min(rxn1_line + 20, len(lines))):
    #     if lines[a].strip() == 'DUPLICATE':
    #         return True
    #     elif lines[a].strip() == '':
    #         return False
    #     elif lines[a][0] != '!' and '=' in lines[a][0]:
    #         return False
    # return False


def mark_as_duplicate(reaction_number):

    assert reaction_number < len(reactions)
    # figure out where to insert

    #print('Marking duplicate:')
    #print('\t', lines[reactions[reaction_number][0]])
    #print()
    rxn1_line = reactions[reaction_number][0]
    for a in range(rxn1_line + 1, min(rxn1_line + 20, len(lines))):
        if lines[a].strip() == '':
            lines.insert(a, 'DUPLICATE\n')
            break
        elif lines[a].strip()[-1] == '/':
            continue
        else:
            lines.insert(a, 'DUPLICATE\n')
            break
    else:
        raise ValueError

    # recalculate reactions and reaction_names
    reaction_number = 0
    for i, line in enumerate(lines):
        if line[0] == '!':
            continue
        trimmed_line = line.split('!')[0]  # remove anything that comes after a comment
        if '=' in trimmed_line:
            reactions[reaction_number] = ((i, line))
            reaction_number += 1


def delete_duplicate_mark(reaction_number):
    # figure out where to delete
    rxn1_line = reactions[reaction_number][0]
    for a in range(rxn1_line, min(rxn1_line + 20, len(lines))):
        if lines[a].strip() == 'DUPLICATE' or 'DUPLICATE' in lines[a].strip().split():
            del(lines[a])
            break

    # recalculate reactions and reaction_names -- happens in 3 places
    reaction_number = 0
    for i, line in enumerate(lines):
        if line[0] == '!':
            continue
        trimmed_line = line.split('!')[0]  # remove anything that comes after a comment
        if '=' in trimmed_line:
            reactions[reaction_number] = ((i, line))
            reaction_number += 1


def comment_out_reaction(reaction_number):

    # if marked_duplicate(reaction_number):
    #     delete_duplicate_mark(reaction_number)
    rxn1_line = reactions[reaction_number][0]
    print(f'modifying: {lines[rxn1_line]}')
    lines[rxn1_line] = '! ' + lines[rxn1_line]

    if 'DUPLICATE' in lines[rxn1_line + 1]:
        lines[rxn1_line + 1] = '! ' + lines[rxn1_line + 1]



# Add missing DUPLICATE markers
for pair in duplicates:
    if not marked_duplicate(pair[0]):

        # delete pair 1
        comment_out_reaction(pair[1])
        print(pair[1])

        # mark_as_duplicate(pair[0])
    # if not marked_duplicate(pair[1]):
    #     mark_as_duplicate(pair[1])

# Get rid of incorrectly marked duplicates
for i in range(0, len(reactions)):
    if marked_duplicate(i) and reaction_names[i] not in duplicate_names:
        delete_duplicate_mark(i)

# need to delete blocks that all start with ! -- only start after the reactions block begins
rxn1_line = reactions[0][0]
del_lines = []
for i in range(rxn1_line, len(lines) - 1):
# for i in range(rxn1_line, 1420):
    
    # find blank line with stuff after it
    if lines[i].strip() == '' and lines[i + 1].strip() != '':
        # print(f'blank {i}  {lines[i]}')
        # print(f'block {i + 1}  {lines[i + 1]}')

        # find the next blank line, checking if all things start with !
        all_comments = True
        for j in range(i + 1, len(lines)):
            if lines[j].strip() == '':
                break
            if lines[j].strip()[0] != '!':
                all_comments = False
        


        if all_comments:
            print(f'Found all comments block on line {i + 1}')
            for k in range(i + 1, j):
                del_lines.append(k)
    
for k in del_lines:
    lines[k] = '\n'
            


    #     for j in range(i + 2, len(lines)):
    #         # find the next blank line after some not blank lines
    #         if lines[j].strip() == '' and lines[j - 1].strip() != '':
    #             # if all the lines are comments, mark them for deletion
    #             # delete_block = True
    #             for k in range(i + 1, j):
    #                 if lines[k].strip() == '':
    #                     print(f'I did this wrong:')
    #                     print(f'i{i}:{lines[i]}')
    #                     print(f'k{k}:{lines[k]}')
    #                     print(f'j{j}:{lines[j]}')
    #                 if lines[k].strip()[0] != '!':
    #                     # delete_block = False
    #                     break
    #             else:
    #                 # delete the block
    #                 for k in range(i + 1, j):
    #                     # lines[k] = '\n'
    #                     del_lines.append(k)
    # for k in del_lines:
    #     lines[k] = '\n'

# save the fixed file
chemkin_file_out = chemkin_file[:-4] + '_fixed' + chemkin_file[-4:]
with open(chemkin_file_out, 'w') as f:
    f.writelines(lines)
