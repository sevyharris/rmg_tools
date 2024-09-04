#script to check how bad my shell guesses are:
import os
import glob
import garbage



i = 915
i = 1714
files = glob.glob(f'/work/westgroup/harris.se/autoscience/autoscience/butane/dft/kinetics/reaction_{i:04}/shell/fwd_ts_*.com')

for j, f in enumerate(files):
    atoms = garbage.read_atoms(f)
    gs = garbage.calculate_garbage_score(atoms)
    #gs = garbage.calculate_garbage_score(atoms, verbose=True)
    print(j, '\t', '{:.3f}'.format(gs))

