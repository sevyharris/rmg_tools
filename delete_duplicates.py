import subprocess
import time
import re

chemkin = 'chem_annotated_fixed.inp'
cmd = f"ck2yaml --input={chemkin} --transport=tran.dat"
cmd_pieces = cmd.split()

b = 0
while True:
    print(f'Attempt {b}')
    proc = subprocess.Popen(cmd_pieces, stdin=None, stdout=subprocess.PIPE, stderr=None, close_fds=True)

    lines = []
    a = 0
    for line in iter(proc.stdout.readline, ''):
        text = line.strip().decode('utf-8')
        if text.lower() == 'passed':
            print('ck2yaml passed!')
            exit(0)
        if text != '':
            a += 1
            lines.append(text)
            print(text)
        if a > 5:
            break

    # scrape the line to delete
    for line in lines:
        if 'found on' in line:
            break
    else:
        raise ValueError('No line found')

    print('extracted this line from ck2yaml:')
    print(line)

    # line = 'found on lines 6849 and 6857 of the kinetics input file.'
    m1=re.search('found on lines (\d*) and (\d*) of the kinetics input file.', line)

    first_line = int(m1.groups()[0])
    second_line = int(m1.groups()[1])
    print(first_line, second_line)


    with open('chem_annotated_fixed.inp', 'r') as f:
        lines = f.readlines()
        print(lines[first_line - 1])
        print(lines[second_line - 1])

    # cound backwards until a space to get the whole comment block
    for i in range(second_line - 1, 0, -1):
        if lines[i].strip() == '':
            break
    # cound forwards until a space to get the whole comment block
    for k in range(second_line, len(lines)):
        if lines[k].strip() == '':
            break
    print()
    print()
    print('This is the comment block:')
    # print the resulting block
    for j in range(i + 1, k):
        print(lines[j])

    # set the resulting block to be empty
    print(f'Deleting the comment block for reaction on line {second_line}')
    for j in range(i + 1, k):
        lines[j] = ''
        
    with open('chem_annotated_fixed.inp', 'w') as f:
        f.writelines(lines)

    print('waiting 2 seconds for next attempt')
    time.sleep(2.0)
    b += 1  
