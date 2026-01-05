# module for fitting kinetics using both fwd and reverse rates
import numpy as np
import rmgpy.kinetics
import scipy.optimize
import copy

def get_reverse_reaction(reaction):
    assert reaction.kinetics is not None
    rev_reaction = copy.deepcopy(reaction)
    tmp_reactants = rev_reaction.reactants
    rev_reaction.reactants = rev_reaction.products
    rev_reaction.products = tmp_reactants
    try:
        rev_reaction.kinetics = reaction.generate_reverse_rate_coefficient()
    except ValueError:
        rev_reaction.kinetics = reaction.generate_reverse_rate_coefficient(surface_site_density=2.72e-5)
    return rev_reaction


def get_fitting(input_rxn, Ts, log10ks, log10ks_rev, fwd_weight=0.5):
    # uses thermo from input_rxn but not kinetics
    # input_rxn is an RMG reaction modified in place and returned with the fitted kinetics
    # Ts are the temperatures at which the rate coefficients are evaluated
    # log10ks are the log10 of the forward rate coefficients at those temperatures
    # log10ks_rev are the log10 of the reverse rate coefficients at those temperatures
    # fwd_weight is the weight given to the forward rate coefficients in the cost function (between 0 and 1)
    
    def cost_function(x):
        log10A, b, Ea = x
        my_rxn_copy = copy.deepcopy(input_rxn)
        test_log10ks = np.zeros_like(Ts)
        rev_test_log10ks = np.zeros_like(Ts)
        my_rxn_copy.kinetics = rmgpy.kinetics.SurfaceArrhenius(A=(np.float_power(10.0, log10A), 'cm^2/(mol*s)'), n=b, Ea=(Ea, 'cal/mol'))
        for i in range(len(Ts)):
            test_log10ks[i] = np.log10(my_rxn_copy.get_rate_coefficient(Ts[i], 101325, surface_site_density=2.4830E-05))
            rev_test_log10ks[i] = np.log10(my_rxn_copy.generate_reverse_rate_coefficient(surface_site_density=2.4830E-05).get_rate_coefficient(Ts[i], 101325))
        error = fwd_weight * np.sum(np.float_power(test_log10ks - log10ks, 2.0)) + (1.0 - fwd_weight) * np.sum(np.float_power(rev_test_log10ks - log10ks_rev, 2.0))
        return error
    
    sol = scipy.optimize.minimize(cost_function, [100, 0.0, 10000], method='Powell',
                             bounds=[(None, None),(None, None),(0.0, None)])  # force positive activation energy

    log10A, b, Ea = sol.x
    fitted_kinetics = rmgpy.kinetics.SurfaceArrhenius(A=(np.float_power(10.0, log10A), 'cm^2/(mol*s)'), n=b, Ea=(Ea, 'cal/mol'))
    input_rxn.kinetics = fitted_kinetics
    return input_rxn