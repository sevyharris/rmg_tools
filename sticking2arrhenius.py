import copy
import rmgpy.kinetics


SDEN = 2.72e-9  # mol/m^2
def sticking2surface_arrhenius(stick_rxn):
    new_rxn = copy.deepcopy(stick_rxn)
    Tlist = np.linspace(300, 3000, 1000)
    klist = np.zeros(len(Tlist))
    for i, T in enumerate(Tlist):
        klist[i] = stick_rxn.get_rate_coefficient(T, 101325, surface_site_density=SDEN) #Value is returned in combination of [m,mol,s]
    new_rxn.kinetics = rmgpy.kinetics.SurfaceArrhenius((1, 'm^5 / mol^2 / s'), 0, (0, 'kcal/mol'))  # dummy reaction
    new_rxn.kinetics.fit_to_data(Tlist, klist, 'm^5 / mol^2 / s')

    return new_rxn
