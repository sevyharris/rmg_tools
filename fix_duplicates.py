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
    if line[0] != '!' and '=' in line:
        reactions.append((i, line))

# print(reactions)
reaction_names = [rxn[1].split()[0] for rxn in reactions]

duplicate_names = set()
duplicates = []
for i in range(0, len(reaction_names)):
    for j in range(0, i):
        if reaction_names[i] == reaction_names[j]:
            duplicates.append((j, i))
            duplicate_names.add(reaction_names[j])


def marked_duplicate(reaction_number):
    rxn1_line = reactions[reaction_number][0]
    for a in range(rxn1_line, min(rxn1_line + 20, len(lines))):
        if lines[a].strip() == 'DUPLICATE':
            return True
        elif lines[a].strip() == '':
            return False
        elif lines[a][0] != '!' and '=' in lines[a][0]:
            return False
    return False


def mark_as_duplicate(reaction_number):
    # figure out where to insert
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

    # recalculate reactions and reaction_names
    reaction_number = 0
    for i, line in enumerate(lines):
        if line[0] != '!' and '=' in line:
            reactions[reaction_number] = ((i, line))
            reaction_number += 1


def delete_duplicate_mark(reaction_number):
    # figure out where to delete
    rxn1_line = reactions[reaction_number][0]
    for a in range(rxn1_line, min(rxn1_line + 20, len(lines))):
        if lines[a].strip() == 'DUPLICATE' or 'DUPLICATE' in lines[a].strip().split():
            del(lines[a])
            break

    # recalculate reactions and reaction_names
    reaction_number = 0
    for i, line in enumerate(lines):
        if line[0] != '!' and '=' in line:
            reactions[reaction_number] = ((i, line))
            reaction_number += 1


# Add missing DUPLICATE markers
for pair in duplicates:
    if not marked_duplicate(pair[0]):
        mark_as_duplicate(pair[0])
    if not marked_duplicate(pair[1]):
        mark_as_duplicate(pair[1])

# Get rid of incorrectly marked duplicates
for i in range(0, len(reactions)):
    if marked_duplicate(i) and reaction_names[i] not in duplicate_names:
        delete_duplicate_mark(i)

# save the fixed file
chemkin_file_out = chemkin_file[:-4] + '_fixed' + chemkin_file[-4:]
with open(chemkin_file_out, 'w') as f:
    f.writelines(lines)
