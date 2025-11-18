# a module to run fast approximations of ignition delay
# using RMG's SimpleReactor
# This is not a real simulation of ignition delay because we're assuming constant T and P
# but perhaps it can be used to get rough estimates quickly


import rmgpy.chemkin
import os
import numpy as np
import matplotlib.pyplot as plt
import rmgpy.solver.simple

import rmgpy.quantity
import rmgpy.rmg.settings


def get_i_thing(thing, thing_list):
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            return i
    return -1


def run_simulation(
    # input_chemkin: str,
    # input_species_dict: str,
    species_list: list,
    reaction_list: list,
    T: float,
    P: float,
    # X: dict | str,
    total_time: float = 1.0,
):
    # # Load the mechanism
    # chemkin = input_chemkin
    # species_dict = input_species_dict

    # species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin, species_dict)

    i_butane = get_i_thing(rmgpy.species.Species(smiles='CCCC'), species_list)
    assert i_butane >= 0
    i_O2 = get_i_thing(rmgpy.species.Species(smiles='[O][O]'), species_list)
    assert i_O2 >= 0
    i_Ar = get_i_thing(rmgpy.species.Species(smiles='[Ar]'), species_list)
    assert i_Ar >= 0


    T = rmgpy.quantity.Quantity(T)
    P = rmgpy.quantity.Quantity(P)
    initial_mole_fractions = {
        species_list[i_butane]: 0.03135,
        species_list[i_Ar]: 0.7649,
        species_list[i_O2]: 0.2038,
    }
    termination_time = rmgpy.solver.TerminationTime(rmgpy.quantity.Quantity(1e1))

    simple_reactor = rmgpy.solver.SimpleReactor(
        T=T,
        P=P,
        initial_mole_fractions=initial_mole_fractions,
        n_sims=1,
        termination=[termination_time],
    )

    # Run the simulation
    model_settings = rmgpy.rmg.settings.ModelSettings(tol_keep_in_edge=0, tol_move_to_core=1, tol_interrupt_simulation=1)
    simulator_settings = rmgpy.rmg.settings.SimulatorSettings()
    terminated, _, invalid_objects, surface_species, surface_reactions, end_time, conversion = simple_reactor.simulate(
        core_species=species_list,
        core_reactions=reaction_list,
        edge_species=[],
        edge_reactions=[],
        surface_species=[],
        surface_reactions=[],
        model_settings=model_settings,
        simulator_settings=simulator_settings,
    )

    # return the simple reactor object
    return simple_reactor


def get_ignition_delay(simple_reactor: rmgpy.solver.SimpleReactor, i_OH: int) -> float:
    # get the ignition delay time from the simple reactor object

    times = [snapshot[0] for snapshot in simple_reactor.snapshots]
    # volumes = [snapshot[1] for snapshot in simple_reactor.snapshots]
    OHs = [snapshot[2 + i_OH] for snapshot in simple_reactor.snapshots]

    # Fine the time at which OH concentration is maximum
    max_OH_index = np.argmax(OHs)
    max_OH_time = times[max_OH_index]
    ignition_delay = max_OH_time

    # find the time at which OH slope is maximum AND dt is greater than 1e-9
    dt = np.diff(times)
    # get first index where dt > 1e-12
    valid_index = np.where(dt > 1e-12)[0]
    if len(valid_index) == 0:
        return ignition_delay  # fallback to max OH time if no valid dt found
    # print(times[valid_index])
    valid_index = int(valid_index[0])
    valid_times = np.array(times[valid_index:])
    valid_OHs = np.array(OHs[valid_index:])

    # slopes = np.gradient(valid_OHs, valid_times)
    # delay_i = np.argmax(slopes)
    # ignition_delay = valid_times[delay_i]

    slopes = np.gradient(valid_OHs, np.log(valid_times))
    delay_i = np.argmax(slopes)
    ignition_delay = valid_times[delay_i]


    # slopes = np.gradient(OHs, times)
    # delay_i = np.argmax(slopes)
    # ignition_delay = times[delay_i]

    return ignition_delay
