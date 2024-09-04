import os
import sys
import glob
import time
import cantera as ct
import numpy as np
import pandas as pd
import concurrent.futures
import rmgpy.chemkin
import subprocess
import scipy
import scipy.optimize
from matplotlib import pyplot as plt
import matplotlib


# TODO - get reaction start index and chunk size to help parallelize...

if len(sys.argv) > 1:
    cti_path = sys.argv[1].replace('.inp', '.yaml')
else:
    raise ValueError("Must provide path to cti file")
aramco = False

working_dir = os.path.join(os.path.dirname(cti_path))
flame_speed_dir = os.path.join(working_dir, 'flame_speeds')
os.makedirs(flame_speed_dir, exist_ok=True)
analyze_convergence = True


def extrapolate_uncertainty(grids, speeds, plot=False):
    """
    Given a list of grid sizes and a corresponding list of flame speeds,
    extrapolate and estimate the uncertainty in the final flame speed.
    Also makes a plot, unless called with `plot=False`.
    """
    grids = list(grids)
    speeds = list(speeds)

    def speed_from_grid_size(grid_size, true_speed, error):
        """
        Given a grid size (or an array or list of grid sizes)
        return a prediction (or array of predictions)
        of the computed flame speed, based on
        the parameters `true_speed` and `error`.

        It seems, from experience, that error scales roughly with
        1/grid_size, so we assume that form.
        """
        return true_speed + error / np.array(grid_size)

    # Fit the chosen form of speed_from_grid_size, to the last four
    # speed and grid size values.
    popt, pcov = scipy.optimize.curve_fit(speed_from_grid_size, grids[-4:], speeds[-4:])

    # How bad the fit was gives you some error, `percent_error_in_true_speed`.
    perr = np.sqrt(np.diag(pcov))
    true_speed_estimate = popt[0]
    percent_error_in_true_speed = perr[0] / popt[0]
    print(
        f"Fitted true_speed is {popt[0] * 100:.4f} ± {perr[0] * 100:.4f} cm/s "
        f"({percent_error_in_true_speed:.1%})"
    )

    # How far your extrapolated infinite grid value is from your extrapolated
    # (or interpolated) final grid value, gives you some other error, `estimated_percent_error`
    estimated_percent_error = (
        speed_from_grid_size(grids[-1], *popt) - true_speed_estimate
    ) / true_speed_estimate
    print(f"Estimated error in final calculation {estimated_percent_error:.1%}")

    # The total estimated error is the sum of these two errors.
    total_percent_error_estimate = abs(percent_error_in_true_speed) + abs(
        estimated_percent_error
    )
    print(f"Estimated total error {total_percent_error_estimate:.1%}")

    if plot:
        plt.clf()
        plt.semilogx(grids, speeds, "o-")
        plt.ylim(
            min(speeds[-5:] + [true_speed_estimate - perr[0]]) * 0.95,
            max(speeds[-5:] + [true_speed_estimate + perr[0]]) * 1.05,
        )
        plt.plot(grids[-4:], speeds[-4:], "or")
        extrapolated_grids = grids + [grids[-1] * i for i in range(2, 8)]
        plt.plot(
            extrapolated_grids, speed_from_grid_size(extrapolated_grids, *popt), ":r"
        )
        plt.xlim(*plt.xlim())
        plt.hlines(true_speed_estimate, *plt.xlim(), colors="r", linestyles="dashed")
        plt.hlines(
            true_speed_estimate + perr[0],
            *plt.xlim(),
            colors="r",
            linestyles="dashed",
            alpha=0.3,
        )
        plt.hlines(
            true_speed_estimate - perr[0],
            *plt.xlim(),
            colors="r",
            linestyles="dashed",
            alpha=0.3,
        )
        plt.fill_between(
            plt.xlim(),
            true_speed_estimate - perr[0],
            true_speed_estimate + perr[0],
            facecolor="red",
            alpha=0.1,
        )

        above = popt[1] / abs(
            popt[1]
        )  # will be +1 if approach from above or -1 if approach from below

        plt.annotate(
            "",
            xy=(grids[-1], true_speed_estimate),
            xycoords="data",
            xytext=(grids[-1], speed_from_grid_size(grids[-1], *popt)),
            textcoords="data",
            arrowprops=dict(
                arrowstyle="|-|, widthA=0.5, widthB=0.5",
                linewidth=1,
                connectionstyle="arc3",
                color="black",
                shrinkA=0,
                shrinkB=0,
            ),
        )

        plt.annotate(
            f"{abs(estimated_percent_error):.1%}",
            xy=(grids[-1], speed_from_grid_size(grids[-1], *popt)),
            xycoords="data",
            xytext=(5, 15 * above),
            va="center",
            textcoords="offset points",
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
        )

        plt.annotate(
            "",
            xy=(grids[-1] * 4, true_speed_estimate - (above * perr[0])),
            xycoords="data",
            xytext=(grids[-1] * 4, true_speed_estimate),
            textcoords="data",
            arrowprops=dict(
                arrowstyle="|-|, widthA=0.5, widthB=0.5",
                linewidth=1,
                connectionstyle="arc3",
                color="black",
                shrinkA=0,
                shrinkB=0,
            ),
        )
        plt.annotate(
            f"{abs(percent_error_in_true_speed):.1%}",
            xy=(grids[-1] * 4, true_speed_estimate - (above * perr[0])),
            xycoords="data",
            xytext=(5, -15 * above),
            va="center",
            textcoords="offset points",
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
        )

        plt.ylabel("Flame speed (m/s)")
        plt.xlabel("Grid size")
        # plt.show()
        plt.savefig(f'convergence_plot_{grid_size}.png')

    return true_speed_estimate, total_percent_error_estimate


def make_callback(flame):
    """
    Create and return a callback function that you will attach to
    a flame solver. The reason we define a function to make the callback function,
    instead of just defining the callback function, is so that it can store
    a pair of lists that persist between function calls, to store the
    values of grid size and flame speed.

    This factory returns the callback function, and the two lists:
    (callback, speeds, grids)
    """
    speeds = []
    grids = []

    def callback(_):
        speed = flame.velocity[0]
        grid = len(flame.grid)
        speeds.append(speed)
        grids.append(grid)
        print(f"Iteration {len(grids)}")
        print(f"Current flame speed is is {speed * 100:.4f} cm/s")
        if len(grids) < 5:
            return 1.0  #
        try:
            extrapolate_uncertainty(grids, speeds)
        except Exception as e:
            print("Couldn't estimate uncertainty. " + str(e))
            return 1.0  # continue anyway
        return 1.0

    return callback, speeds, grids


# perturb every species and reaction in the mechanism
# we'll select the perturbations one at a time later in the script
def perturb_species(species, delta):
    # takes in an RMG species object
    # change the enthalpy offset
    increase = None
    for poly in species.thermo.polynomials:
        new_coeffs = poly.coeffs
        if not increase:
            # Only define the increase in enthalpy once or you'll end up with numerical gaps in continuity
            increase = delta * new_coeffs[5]
        new_coeffs[5] += increase
        poly.coeffs = new_coeffs


def perturb_reaction(rxn, delta):  # 0.1 is default
    # takes in an RMG reaction object
    # delta is the ln(k) amount to perturb the A factor
    # delta is a multiplicative factor- units don't matter, yay!
    # does not deepycopy because there's some issues with rmgpy.reactions copying
    if type(rxn.kinetics) == rmgpy.kinetics.chebyshev.Chebyshev:
        rxn.kinetics.coeffs.value_si[0][0] += np.log10(1.0 + delta)
    elif type(rxn.kinetics) in [rmgpy.kinetics.falloff.Troe, rmgpy.kinetics.falloff.ThirdBody, rmgpy.kinetics.falloff.Lindemann]:
        if hasattr(rxn.kinetics, 'arrheniusHigh'):
            rxn.kinetics.arrheniusHigh.A.value *= np.exp(delta)
        if hasattr(rxn.kinetics, 'arrheniusLow'):
            rxn.kinetics.arrheniusLow.A.value *= np.exp(delta)
    elif type(rxn.kinetics) == rmgpy.kinetics.arrhenius.MultiArrhenius:
        for j in range(len(rxn.kinetics.arrhenius)):
            rxn.kinetics.arrhenius[j].A.value *= np.exp(delta)
    elif type(rxn.kinetics) == rmgpy.kinetics.arrhenius.PDepArrhenius:
        for j in range(len(rxn.kinetics.arrhenius)):
            if type(rxn.kinetics.arrhenius[j]) == rmgpy.kinetics.arrhenius.Arrhenius:
                rxn.kinetics.arrhenius[j].A.value *= np.exp(delta)
            elif type(rxn.kinetics.arrhenius[j]) == rmgpy.kinetics.arrhenius.MultiArrhenius:
                for k in range(len(rxn.kinetics.arrhenius[j].arrhenius)):
                    rxn.kinetics.arrhenius[j].arrhenius[k].A.value *= np.exp(delta)
            else:
                raise ValueError(f'weird kinetics {str(rxn.kinetics)}')
    elif type(rxn.kinetics) == rmgpy.kinetics.arrhenius.MultiPDepArrhenius:
        for i in range(len(rxn.kinetics.arrhenius)):
            for j in range(len(rxn.kinetics.arrhenius[i].arrhenius)):
                if type(rxn.kinetics.arrhenius[i].arrhenius[j]) == rmgpy.kinetics.arrhenius.Arrhenius:
                    rxn.kinetics.arrhenius[i].arrhenius[j].A.value *= np.exp(delta)
                elif type(rxn.kinetics.arrhenius[i].arrhenius[j]) == rmgpy.kinetics.arrhenius.MultiArrhenius:
                    for k in range(len(rxn.kinetics.arrhenius[i].arrhenius[j].arrhenius)):
                        rxn.kinetics.arrhenius[i].arrhenius[j].arrhenius[k].A.value *= np.exp(delta)
                else:
                    raise ValueError(f'weird kinetics {str(rxn.kinetics)}')

    else:  # Arrhenius
        rxn.kinetics.A.value *= np.exp(delta)


base_yaml_path = os.path.join(working_dir, 'base.yaml')
perturbed_chemkin = os.path.join(working_dir, 'perturbed.inp')
perturbed_yaml_path = os.path.join(working_dir, 'perturbed.yaml')

skip_create_perturb = False
if os.path.exists(perturbed_yaml_path):
    skip_create_perturb = True
    print('Perturbed yaml already exists, skipping creation of perturbed mechanism')

if not skip_create_perturb:
    # load the chemkin file and create a normal and perturbed cti for simulations:
    # # write base cantera
    chemkin = cti_path.replace('.yaml', '.inp')
    transport = os.path.join(working_dir, 'tran.dat')
    species_dict = os.path.join(working_dir, 'species_dictionary.txt')

    species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(chemkin, dictionary_path=species_dict, transport_path=transport, use_chemkin_names=True)
    print(f'Loaded {len(species_list)} species, {len(reaction_list)} reactions')

    subprocess.run(['ck2yaml', f'--input={chemkin}', f'--transport={transport}', f'--output={base_yaml_path}'])

    delta = 0.1
    for i in range(0, len(species_list)):
        perturb_species(species_list[i], delta)

    for i in range(0, len(reaction_list)):
        try:
            perturb_reaction(reaction_list[i], delta)
        except AttributeError:
            continue

    # save the results
    rmgpy.chemkin.save_chemkin_file(perturbed_chemkin, species_list, reaction_list, verbose=True, check_for_duplicates=True)
    subprocess.run(['ck2yaml', f'--input={perturbed_chemkin}', f'--transport={transport}', f'--output={perturbed_yaml_path}'])


# load the 2 ctis
base_gas = ct.Solution(base_yaml_path)
perturbed_gas = ct.Solution(perturbed_yaml_path)


# load the experimental conditions
flame_speed_data = '/work/westgroup/harris.se/autoscience/autoscience/butane/experimental_data/butane_flamespeeds.csv'
# flame_speed_data = '/home/moon/autoscience/autoscience/butane/experimental_data/butane_flamespeeds.csv'
df_exp = pd.read_csv(flame_speed_data)

initial_guess_path = os.path.dirname(cti_path)

# get just the Park data
data_slice = df_exp[df_exp['Reference'] == 'Park et al. 2016']
N = 51
phi_min = 0.6
phi_max = 2.0
equiv_ratios = np.linspace(phi_min, phi_max, N)
temperatures = np.ones(len(equiv_ratios)) * data_slice['Tu (K)'].values[0]
pressures = np.ones(len(equiv_ratios)) * data_slice['Pu (atm)'].values[0] * ct.one_atm

# list of starting conditions
# Define stoichiometric coefficients
v_fuel = 1.0
v_oxidizer = 13.0 / 2.0
v_N2 = 0.79 * (v_oxidizer / 0.21)  # air is approximately 79% N2 and 21% O2

# calculate actual ratio of fuel to oxidizer
actual_ratio = equiv_ratios * (v_fuel / v_oxidizer)
# start with 1.0 oxidizer, then normalize
x_O2 = 1.0
x_C4H10 = actual_ratio * x_O2
x_N2 = 0.79 * (x_O2 / .21)
total = x_O2 + x_C4H10 + x_N2
x_O2 = x_O2 / total
x_C4H10 = x_C4H10 / total
x_N2 = x_N2 / total

concentrations = [{'butane(1)': x_C4H10[i], 'O2(2)': x_O2[i], 'N2': x_N2[i]} for i in range(0, len(equiv_ratios))]


def get_max_index_from_globlist(globlist):
    # Function for making use of previous flame runs
    FINDSTR = 'saved_flame_'
    underscore_loc = globlist[0].find(FINDSTR) + len(FINDSTR)
    max_idx = -1
    for g in globlist:
        if int(g[underscore_loc: -3]) > max_idx:
            max_idx = int(g[underscore_loc: -3])
    return max_idx


# function for running a flame speed
# assumes gas has been properly initialized
def run_flame_speed(condition_index):
    gas = ct.Solution(cti_path)

    if 'butane(1)' not in gas.species_names and 'butane(1)' in concentrations[condition_index].keys():
        # check the keys to make sure we don't pop twice TODO should be O2 not in keys
        concentrations[condition_index]['C4H10'] = concentrations[condition_index].pop('butane(1)')
    if 'O2(2)' not in gas.species_names and 'O2(2)' in concentrations[condition_index].keys():
        concentrations[condition_index]['O2'] = concentrations[condition_index].pop('O2(2)')
    if 'CO2(7)' not in gas.species_names and 'CO2(7)' in concentrations[condition_index].keys():
        concentrations[condition_index]['CO2'] = concentrations[condition_index].pop('CO2(7)')
    if 'Ar' not in gas.species_names and 'Ar' in concentrations[condition_index].keys():
        concentrations[condition_index]['AR'] = concentrations[condition_index].pop('Ar')

    gas.TPX = temperatures[condition_index], pressures[condition_index], concentrations[condition_index]

    tol_ss = [1.0e-13, 1.0e-9]  # abs and rel tolerances for steady state problem
    tol_ts = [1.0e-13, 1.0e-9]  # abs and rel tie tolerances for time step function

    width = 0.08
    flame = ct.FreeFlame(gas, width=width)
    flame.flame.set_steady_tolerances(default=tol_ss)   # set tolerances
    flame.flame.set_transient_tolerances(default=tol_ts)
    flame.set_refine_criteria(ratio=5, slope=0.25, curve=0.27)  # decrease this
    # flame.max_time_step_count = 5000
    flame.max_time_step_count = 900
    loglevel = 1

    # set up from previous run
    h5_filepath_guess = os.path.join(flame_speed_dir, f"saved_rxn_flame_{condition_index}.h5")
    if not os.path.exists(h5_filepath_guess):
        # Attempt to reuse base solves
        # MAKE SURE THE saved_flame_{} has the same index. Need to convert to closest thing.
        globlist = glob.glob(os.path.join(flame_speed_dir, 'saved_flame_*.h5'))
        base_N = get_max_index_from_globlist(globlist)
        N_to_use = int(condition_index / N * base_N)
        h5_filepath_guess = os.path.join(flame_speed_dir, f"saved_flame_{N_to_use}.h5")  # can try a guess from regular version

    if os.path.exists(h5_filepath_guess):

        print('Loading initial flame conditions from', h5_filepath_guess)
        print("Load initial guess from HDF file via SolutionArray")
        arr2 = ct.SolutionArray(gas)
        # the flame domain needs to be specified as subgroup
        arr2.read_hdf(h5_filepath_guess, group="freeflame", subgroup="flame", force=True)
        flame.set_initial_guess(data=arr2)

    else:
        print('Initial flame conditions not found, using default')
        # raise OSError

    if analyze_convergence:
        callback, speeds, grids = make_callback(flame)
        flame.set_steady_callback(callback)

    print("about to solve")
    flame.solve(loglevel=loglevel, auto=True)
    Su = flame.velocity[0]

    print("Save HDF")  # Maybe a bad idea to save when you're reusing across species... saving could lead to a HUGE number of files
    hdf_filepath = os.path.join(flame_speed_dir, f"saved_rxn_flame_{condition_index}.h5")
    flame.write_hdf(
        hdf_filepath,
        group="freeflame",
        mode="w",
        quiet=False,
        description=("butane flame"),
    )

    return Su


# compute and save the delays
reaction_flame_speeds = np.zeros((len(perturbed_gas.reactions()), len(equiv_ratios)))

for i in range(0, len(perturbed_gas.reactions())):
    print(f'perturbing {i} {perturbed_gas.reactions()[i]}')
    # load the base gas
    base_gas = ct.Solution(base_yaml_path)

    if same_reaction(perturbed_gas.reactions()[i], base_gas.reactions()[i]):
        perturbed_index = i
    else:
        for k in range(0, len(base_gas.reactions())):
            if same_reaction(perturbed_gas.reactions()[k], base_gas.reactions()[i]):
                perturbed_index = k
                break
        else:
            raise ValueError('Could not find matching reaction in base mechanism')

    base_gas.modify_reaction(i, perturbed_gas.reactions()[perturbed_index])

    # Run all simulations in parallel
    flame_speeds = np.zeros(len(equiv_ratios))
    condition_indices = np.arange(0, len(equiv_ratios))

    with concurrent.futures.ProcessPoolExecutor(max_workers=16) as executor:
        for condition_index, flame_speed in zip(condition_indices, executor.map(run_flame_speed, condition_indices)):
            flame_speeds[condition_index] = flame_speed

    reaction_flame_speeds[i, :] = flame_speeds

# save the result as a pandas dataframe
np.save(os.path.join(flame_speed_dir, f'reaction_flame_speeds.npy'), reaction_flame_speeds)
