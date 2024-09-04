# script to sanity check the geometry of a molecule and tell you if it's garbage
import sys
import ase.io
import ase.io.gaussian


def calculate_garbage_score(atoms, verbose=False):
    garbage_score = 0
    min_H_dist = 0.75
    min_heavy_dist = 1.1
    nn_dist = 1.5
    nn_threshold = 4

    # do a nearest neighbor check to make sure no atom is too close to another
    num_combos = len(atoms) * (len(atoms) - 1) / 2
    for i in range(0, len(atoms)):
        # count the atoms within 1.5A of this atom
        n_close = 0

        for j in range(i + 1, len(atoms)):
            dist = atoms.get_distance(i, j)
            if atoms[i].symbol == 'H' or atoms[j].symbol == 'H':
                threshold = min_H_dist
            else:
                threshold = min_heavy_dist
            if dist < threshold:
                atomic_weight = atoms.get_atomic_numbers()[i] + atoms.get_atomic_numbers()[j]
                garbage_score += atomic_weight / num_combos + threshold - dist
                if verbose:
                    print(f'Atom {i} and atom {j} are too close ({dist:.3f} < {threshold:.3f})')
            if dist < nn_dist:
                n_close += 1
        for j in range(0, i):
            dist = atoms.get_distance(i, j)
            if dist < nn_dist:
                n_close += 1

        if n_close > nn_threshold:
            garbage_score += n_close / 10.0
            if verbose:
                print(f'Atom {i} has {n_close} neighbors')
    return garbage_score


def read_atoms(input_file):
    # convert the geometry into ase form
    if input_file.endswith('.log'):
        with open(input_file, 'r') as f:
            atoms = ase.io.gaussian.read_gaussian_out(f)

    elif input_file.endswith('.com') or input_file.endswith('.gjf'):
        try:
            with open(input_file, 'r') as f:
                atoms = ase.io.gaussian.read_gaussian_in(f)
        except ase.io.ParseError:
            with open(input_file, 'r') as f:
                # code copied from ase.io.gaussian
                # https://wiki.fysik.dtu.dk/ase/_modules/ase/io/gaussian.html

                parameters = {}
                file_sections = ase.io.gaussian._get_gaussian_in_sections(f)
                parameters.update(ase.io.gaussian._get_all_link0_params(file_sections['link0']))
                parameters.update(ase.io.gaussian._get_all_route_params(file_sections['route']))
                parameters.update(ase.io.gaussian._get_charge_mult(file_sections['charge_mult']))
                atoms, _ = ase.io.gaussian._get_atoms_from_molspec(
                    file_sections['mol_spec'])
                gc = ase.io.gaussian.GaussianConfiguration(atoms, parameters)
                atoms = gc.get_atoms()
    elif input_file.endswith('.xyz'):
        raise NotImplementedError('xyz files not yet supported')
    else:
        raise ValueError(f'Unrecognized file extension for {input_file}')
    return atoms


if __name__ == '__main__':
    if len(sys.argv) < 2:
        # input_file = '/home/moon/autoscience/reaction_calculator/sanity_check_examples/fwd_ts_0000.log'
        # input_file = '/home/moon/autoscience/reaction_calculator/sanity_check_examples/fwd_ts_0000.com'
        # input_file = '/home/moon/autoscience/reaction_calculator/sanity_check_examples/fwd_ts_0019.log'
        input_file = '/home/moon/autoscience/reaction_calculator/sanity_check_examples/fwd_ts_0019.com'
    else:
        input_file = sys.argv[1]
    atoms = read_atoms(input_file)
    garbage_score = calculate_garbage_score(atoms)
    print(f'Garbage score: {garbage_score:.3f}')
