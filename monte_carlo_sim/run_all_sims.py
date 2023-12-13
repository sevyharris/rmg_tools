# script to set up many simulations
import random
# import pandas as pd
# import pickle
import numpy as np
import cantera as ct
import concurrent.futures

import rmgpy.chemkin
# from rmgpy.species import Species
# from rmgpy.tools.canteramodel import Cantera, get_rmg_species_from_user_species
# # from rmgpy.tools.globaluncertainty import ReactorPCEFactory
# from rmgpy.tools.globaldelayuncertainty import ReactorPCEFactory
# import rmgpy.tools.uncertainty_gao


# decide which chunk of the random sampling to handle
sample_start_index = 0


# first, read in the model
chemkin_file = '/home/moon/hw8/bill_green_kinetics/project/base24_analysis/chem_annotated.inp'
dict_file = '/home/moon/hw8/bill_green_kinetics/project/base24_analysis/species_dictionary.txt'
transport = '/home/moon/hw8/bill_green_kinetics/project/base24_analysis/tran.dat'

cti_path = '/home/moon/hw8/bill_green_kinetics/project/base24_analysis/chem_annotated.yaml'


# load the uncertainty values so we know how much to perturb each parameter
species_uncertainty_file = '/home/moon/hw8/bill_green_kinetics/project/base24_analysis/delta_Gs.npy'
reaction_uncertainty_file = '/home/moon/hw8/bill_green_kinetics/project/base24_analysis/delta_ks.npy'
delta_ks = np.load(reaction_uncertainty_file)
delta_Gs = np.load(species_uncertainty_file)


species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin_file, dictionary_path=dict_file, transport_path=transport)


# experimental conditions copied from subset of # https://doi-org.ezproxy.neu.edu/10.1016/j.combustflame.2010.01.016
conc_dicts = [{
    'O2(2)': 0.2038,
    'butane(1)': 0.03135,
    'Ar': 0.7649
}] * 10
experimental_Ts = [700, 731, 756, 794, 840, 872, 906, 926, 959, 1000]
experimental_Ps = [1013250] * 10



# generate the perturbation matrix...


# Take Reactor Conditions from Table 7 of supplementary info in
# https://doi-org.ezproxy.neu.edu/10.1016/j.combustflame.2010.01.016
def get_delay(gas, T, P, X):
    # function to run a RCM simulation

    t_end = 1.0  # time in seconds
    gas.TPX = T, P, X

    env = ct.Reservoir(ct.Solution('air.yaml'))
    # env = ct.Reservoir(ct.Solution('air.xml'))
    reactor = ct.IdealGasReactor(gas)
    wall = ct.Wall(reactor, env, A=1.0, velocity=0)
    reactor_net = ct.ReactorNet([reactor])
    # # allegedly faster solving
    # reactor_net.derivative_settings = {"skip-third-bodies": True, "skip-falloff": True}
    # reactor_net.preconditioner = ct.AdaptivePreconditioner()

    times = [0]
    T = [reactor.T]
    P = [reactor.thermo.P]
    X = [reactor.thermo.X]  # mol fractions
    while reactor_net.time < t_end:
        reactor_net.step()

        times.append(reactor_net.time)
        T.append(reactor.T)
        P.append(reactor.thermo.P)
        X.append(reactor.thermo.X)

    slopes = np.gradient(P, times)
    i = np.argmax(slopes)
    return times[i]


def run_simulation(condition_index):
    gas = ct.Solution(cti_path)
    X = conc_dicts[condition_index]
    delay = get_delay(gas, experimental_Ts[condition_index], experimental_Ps[condition_index], X)
    # print(f'Completed {condition_index}:\t {delay}')
    return delay

N = 2
total_delays = np.zeros((N, len(experimental_Ts)))
for i in range(0, N):  # do 1000 samples

    # Run 16 simulations in parallel
    delays = np.zeros(len(experimental_Ts))
    condition_indices = np.arange(0, len(experimental_Ts))
    with concurrent.futures.ProcessPoolExecutor(max_workers=len(experimental_Ts)) as executor:
        for condition_index, delay_time in zip(condition_indices, executor.map(run_simulation, condition_indices)):
            delays[condition_index] = delay_time
    total_delays[i, :] = delays

np.save(f'delays_{sample_start_index:04}.npy', total_delays)
