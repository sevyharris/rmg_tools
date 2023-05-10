# script to report the lowest energy TS of the given logs
import sys
import glob
import arkane.ess.gaussian


log_search_string = sys.argv[1]
log_files = glob.glob(log_search_string)

energies = []
valid_TS = []

best_option = None
lowest_energy = 1e10

for logfile in log_files:
    try:
        g_reader = arkane.ess.gaussian.GaussianLog(logfile)
        freq = g_reader.load_negative_frequency()
    except arkane.exceptions.LogError:
        energies.append(lowest_energy)
        valid_TS.append(False)
        continue

    energy = g_reader.load_energy()
    energies.append(energy)

    # get a list of the vibrational frequencies
    if freq < -500:
        valid_TS.append(True)
        if energy < lowest_energy:
            lowest_energy = energy
            best_option = logfile
    else:
        valid_TS.append(False)

# print the results of each TS
for i, logfile in enumerate(log_files):
    print('TS: {} Energy: {} Valid: {}'.format(logfile, energies[i], valid_TS[i]))

print('The lowest energy TS is: {}'.format(best_option))
