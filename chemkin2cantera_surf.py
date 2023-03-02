# script to convert the chem_annotated-gas.inp and chem_annotated-surf.inp files to a cantera input file
import os
import sys

import rmgpy.chemkin
import cantera.ck2cti


# read in the gas and surf files
if len(sys.argv) != 3:
    gas_ck = '/home/moon/rmg/RMG-tests/tests/rmg/catalytic_combustion/chemkin/chem_annotated-gas.inp'
    surf_ck = '/home/moon/rmg/RMG-tests/tests/rmg/catalytic_combustion/chemkin/chem_annotated-surface.inp'
elif sys.argv[1].endswith('-gas.inp'):
    gas_ck = sys.argv[1]
    surf_ck = sys.argv[2]
else:
    gas_ck = sys.argv[2]
    surf_ck = sys.argv[1]

transport = os.path.join(os.path.dirname(gas_ck), 'tran.dat')
spec_dict = os.path.join(os.path.dirname(gas_ck), 'species_dictionary.txt')

gas_species, gas_reactions = rmgpy.chemkin.load_chemkin_file(gas_ck, dictionary_path=spec_dict, transport_path=transport)
surf_species, surf_reactions = rmgpy.chemkin.load_chemkin_file(surf_ck, dictionary_path=spec_dict, transport_path=transport)


# code copied from rmgpy.main
def generate_cantera_files(chemkin_file, **kwargs):
    """
    Convert a chemkin mechanism chem.inp file to a cantera mechanism file chem.cti
    and save it in the cantera directory
    """
    output_dir = os.path.dirname(chemkin_file)
    transport_file = os.path.join(os.path.dirname(chemkin_file), 'tran.dat')
    file_name = os.path.splitext(os.path.basename(chemkin_file))[0] + '.cti'
    out_name = os.path.join(output_dir, 'cantera', file_name)
    if 'surfaceFile' in kwargs:
        out_name = out_name.replace('-gas.', '.')
    cantera_dir = os.path.dirname(out_name)
    try:
        os.makedirs(cantera_dir)
    except OSError:
        if not os.path.isdir(cantera_dir):
            raise
    if os.path.exists(out_name):
        os.remove(out_name)
    parser = cantera.ck2cti.Parser()
    try:
        parser.convertMech(chemkin_file, transportFile=transport_file, outName=out_name, quiet=True, permissive=True,
                            **kwargs)
    except cantera.ck2cti.InputParseError:
        print("Error converting to Cantera format.")
        print("Trying again without transport data file.")
        parser.convertMech(chemkin_file, outName=out_name, quiet=True, permissive=True, **kwargs)


output_dir = os.path.dirname(gas_ck)
generate_cantera_files(
    os.path.join(output_dir, 'chem-gas.inp'),
    surfaceFile=(os.path.join(output_dir,'chem-surface.inp'))
)
generate_cantera_files(
    os.path.join(output_dir, 'chem_annotated-gas.inp'),
    surfaceFile=(os.path.join(output_dir,'chem_annotated-surface.inp'))
)
