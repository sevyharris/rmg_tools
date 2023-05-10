# script to fix Gaussian error where it lists rotational constant as *******

import sys

gaussian_file = sys.argv[1]
# new_file = gaussian_file[:-4] + '_fixed' + gaussian_file[-4:]
new_file = gaussian_file

lines_with_rot_const = []
with open(gaussian_file, 'r') as f:
    lines = f.readlines()
    # Find the line with the rotational constant
    for i, line in enumerate(lines):
        if 'Rotational constants' in line:
            rot_const_line = i
            lines_with_rot_const.append(i)

if len(lines_with_rot_const) < 2:
    print('No alternate rotational constant found')
    sys.exit()

# Find the line with the asterisks
for k, i in enumerate(lines_with_rot_const):
    line = lines[i]
    if '******' in line:
        # tokenize
        tokens = line.split()
        # find the token with the ****
        for j, token in enumerate(tokens):
            if '******' in token:
                # replace the **** with the previous value
                last_index = lines_with_rot_const[k - 1]
                prev_value = lines[last_index].split()[j]
                tokens[j] = prev_value

                # save the new file
                lines[i] = ' ' + ' '.join(tokens) + '\n'
                with open(new_file, 'w') as f:
                    f.writelines(lines)
                print('New file saved as {}'.format(new_file))
                sys.exit()
else:
    print('No asterisks found')
