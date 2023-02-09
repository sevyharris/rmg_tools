# script to print the energy of a gaussian log file
import sys
import arkane.ess.gaussian


logfile = sys.argv[1]

g_reader = arkane.ess.gaussian.GaussianLog(logfile)
energy = g_reader.load_energy()

print(energy)

