# script to run many flame speeds
import os
import cantera as ct
import numpy as np
import pandas as pd
import concurrent.futures
import sys
from matplotlib import pyplot as plt
import matplotlib
import scipy
import scipy.optimize


if len(sys.argv) > 1:
    cti_path = sys.argv[1]
else:
    raise ValueError("Must provide path to cti file")


def extrapolate_uncertainty(grids, speeds, plot=True):
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


if cti_path.endswith('.inp'):
    cti_path = cti_path.replace('.inp', '.yaml')

gas = ct.Solution(cti_path)

working_dir = os.path.join(os.path.dirname(cti_path), 'flame_speeds')
os.makedirs(working_dir, exist_ok=True)
analyze_convergence = True

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

# concentrations = [{'C4H10': x_C4H10[i], 'O2': x_O2[i], 'N2': x_N2[i]} for i in range(0, len(equiv_ratios))]
concentrations = [{'butane(1)': x_C4H10[i], 'O2(2)': x_O2[i], 'N2': x_N2[i]} for i in range(0, len(equiv_ratios))]


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
    h5_filepath_guess = os.path.join(working_dir, f"saved_flame_{condition_index}.h5")
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

    # print("Save CSV")
    # # csv_filepath = os.path.join(os.path.dirname(cti_path), f"flame_{condition_index}.csv")
    # csv_filepath = os.path.join(working_dir, f"saved_flame_{condition_index}.csv")
    # flame.write_csv(csv_filepath)

    print("Save HDF")
    hdf_filepath = os.path.join(working_dir, f"saved_flame_{condition_index}.h5")
    flame.write_hdf(
        hdf_filepath,
        group="freeflame",
        mode="w",
        quiet=False,
        description=("butane flame"),
    )

    # print("Save YAML")
    # yaml_filepath = os.path.join(working_dir, f"saved_flame_{condition_index}.yaml")
    # flame.save(yaml_filepath, name="solution", description=f"butane flame {condition_index}")

    return Su


# Run all simulations in parallel
flame_speeds = np.zeros(len(equiv_ratios))
condition_indices = np.arange(0, len(equiv_ratios))
with concurrent.futures.ProcessPoolExecutor(max_workers=16) as executor:
    for condition_index, flame_speed in zip(condition_indices, executor.map(run_flame_speed, condition_indices)):
        flame_speeds[condition_index] = flame_speed


results_df = pd.DataFrame(columns=['T', 'P', 'flame speed (cm/s)', 'phi', 'X'])
for i in range(0, len(temperatures)):
    results_df = results_df.append({
        'T': temperatures[i],
        'P': pressures[i],
        'flame speed (cm/s)': flame_speeds[i],
        'phi': equiv_ratios[i],
        'X': concentrations[i],
    }, ignore_index=True)
results_df.to_csv(os.path.join(working_dir, 'flame_speed_results.csv'))


out_df = pd.DataFrame(flame_speeds)
# out_df.to_csv('aramco_flame_speeds.csv')

# out_df.to_csv('naive_improved_flame_speeds.csv')
# out_df.to_csv('rmg_changed_flame_speeds.csv')
# out_df.to_csv(f'{cti_path[:-4]}.csv')
out_df.to_csv(f'{os.path.join(working_dir, os.path.basename(cti_path))[:-4]}.csv')
