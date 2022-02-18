import os
from rmgpy.tools import fluxdiagram



model_dir = '/home/moon/rmg/rmg_tools/uncertainty/nheptane'


chemkin_file = os.path.join(model_dir, 'chem_annotated.inp')
dict_file = os.path.join(model_dir, 'species_dictionary.txt')
input_file = os.path.join(model_dir, 'flux_input.py')
output_path = "/home/moon/rmg/rmg_tools/uncertainty/nheptane/mech_viz"

settings = {
    'max_node_count': 100,
    'max_edge_count': 100,
    # 'concentration_tol': 1e-9,
    # 'species_rate_tol': 1e-9
    'max_node_pen_width': 5,
    'max_edge_pen_width': 5,

}
fluxdiagram.create_flux_diagram(input_file, chemkin_file, dict_file, save_path=output_path, settings=settings)
