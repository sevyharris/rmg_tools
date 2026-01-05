# module for generating flux diagrams for PFR in Cantera


import os

import numpy as np
import pydot
import rmgpy.tools.fluxdiagram

import cantera as ct


def get_i_thing(thing, thing_list):
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            return i
    return -1


# Options controlling the individual flux diagram renderings:
program = 'dot'  # The program to use to lay out the nodes and edges
max_node_count = 42  # The maximum number of nodes to show in the diagram
max_edge_count = 42  # The maximum number of edges to show in the diagram
concentration_tol = 1e-6  # The lowest fractional concentration to show (values below this will appear as zero)
species_rate_tol = 1e-6  # The lowest fractional species rate to show (values below this will appear as zero)
max_node_pen_width = 7.0  # The thickness of the border around a node at maximum concentration
max_edge_pen_width = 9.0  # The thickness of the edge at maximum species rate
radius = 1  # The graph radius to plot around a central species
central_reaction_count = None  # The maximum number of reactions to draw from each central species (None draws all)
# If radius > 1, then this is the number of reactions from every species

# Options controlling the ODE simulations:
initial_time = 1e-12  # The time at which to initiate the simulation, in seconds
time_step = 10 ** 0.1  # The multiplicative factor to use between consecutive time points
abs_tol = 1e-16  # The absolute tolerance to use in the ODE simluations
rel_tol = 1e-8  # The relative tolerance to use in the ODE simulations

# Options controlling the generated movie:
video_fps = 6  # The number of frames per second in the generated movie
initial_padding = 5  # The number of seconds to display the initial fluxes at the start of the video
final_padding = 5  # The number of seconds to display the final fluxes at the end of the video

# species_path = None
java = False            # always False
chemkin_output = ''     # this will be generated automatically
central_species_list = None
superimpose = False     # this will always be false, delete it
save_states = False
read_states = False     # fine to keep this always false and delete relevant code below
# diffusion_limited = True
diffusion_limited = True
check_duplicates = True


def make_flux_diagram_pfr(chemkin_gas, distances1, gas_concs, surf_covs, gas_rates, surf_rates, diagram_base_name='pfr_flux'):


    # rmg_input_file = '/home/moon/nitridation/postdoc/repro_rmg_attempts/v3.1.0/min_input.py'  # for conditions <---- Note only one T,P,X in the reactor
    rmg_input_file = '/home/moon/uncertainty_estimator/cpox_pt/cantera_versions/pfr/min_input.py'

    # diagram_base_name = 'pfr_flux'
    os.makedirs(diagram_base_name, exist_ok=True)
    working_dir = os.path.join(os.getcwd(), diagram_base_name)

    chemkin_surf = chemkin_gas.replace('chem_annotated-gas.inp', 'chem_annotated-surface.inp')
    chemkin_dict = chemkin_gas.replace('chem_annotated-gas.inp', 'species_dictionary.txt')
    
    # mech_1_inp = '/home/moon/uncertainty_estimator/beef/cpox_Pt/chem_annotated-gas.inp'
    # mech_1_surf = '/home/moon/uncertainty_estimator/beef/cpox_Pt/chem_annotated-surface.inp'
    # mech_1_dict = '/home/moon/uncertainty_estimator/beef/cpox_Pt/species_dictionary.txt'
    # mech_1_label = 'debug_fast_kinetics'

    # mech_1_label = 'debug_fast_kinetics_R24_0p15'
    mech_1_label = 'fast_kinetics'


    generate_images = True
    species_path = '/home/moon/nitridation/postdoc/repro_rmg_attempts/v3.1.0/species'
    print('Loading RMG job 1...')
    rmg_job1 = rmgpy.tools.fluxdiagram.load_rmg_job(
        rmg_input_file,
        chemkin_gas,
        chemkin_dict,
        surface_path=chemkin_surf,
        generate_images=generate_images,
        check_duplicates=check_duplicates
    )


    if isinstance(distances1, str):
        distances1 = np.load(distances1)
    if isinstance(gas_concs, str):
        gas_concs = np.load(gas_concs)
    if isinstance(surf_covs, str):
        surf_covs = np.load(surf_covs)
    if isinstance(gas_rates, str):
        gas_rates = np.load(gas_rates)
    if isinstance(surf_rates, str):
        surf_rates = np.load(surf_rates)

    if len(distances1) > 1000:
        DOWNSAMPLE_RATE = int(np.ceil(len(distances1) / 1000.0))
        distances1 = distances1[::DOWNSAMPLE_RATE]
        
        gas_concs = gas_concs[::DOWNSAMPLE_RATE, :]
        surf_covs = surf_covs[::DOWNSAMPLE_RATE, :]
        gas_rates = gas_rates[::DOWNSAMPLE_RATE, :]
        surf_rates = surf_rates[::DOWNSAMPLE_RATE, :]


    concentrations1 = np.concatenate((gas_concs, surf_covs), axis=1)
    reaction_rates1 = np.concatenate((gas_rates, surf_rates), axis=1)
    species_list1 = rmg_job1.reaction_model.core.species[:]
    reaction_list1 = rmg_job1.reaction_model.core.reactions[:]
    num_species1 = len(species_list1)

    def get_fake_pairs(reaction):
        pairs = []
        if reaction.is_isomorphic(reaction_list1[410]): # oxygen abstraction
            pairs.append((reaction.reactants[0], reaction.products[0]))
            pairs.append((reaction.reactants[1], reaction.products[1]))
            pairs.append((reaction.reactants[1], reaction.products[0]))
            return pairs
        
        for r in reaction.reactants:
            for p in reaction.products:
                pairs.append((r, p))
        return pairs
    

    # Compute the rates between each pair of species (build up big matrices)
    species_rates1 = np.zeros((len(distances1), num_species1, num_species1), float)
    for index1, reaction1 in enumerate(reaction_list1):
        rate1 = reaction_rates1[:, index1]
        reaction1.generate_pairs()
        
        
    #     if not reaction1.pairs: reaction1.generate_pairs()
        for reactant1, product1 in reaction1.pairs:
    #     for reactant1, product1 in get_fake_pairs(reaction1):
            reactant_index1 = species_list1.index(reactant1)
            product_index1 = species_list1.index(product1)
            species_rates1[:, reactant_index1, product_index1] += rate1
            species_rates1[:, product_index1, reactant_index1] -= rate1
            

    # Determine the maximum concentration for each species and the maximum overall concentration
    max_concentrations1 = np.max(np.abs(concentrations1), axis=0)
    max_concentration1 = np.max(max_concentrations1)

    # Determine the maximum reaction rates
    max_reaction_rates1 = np.max(np.abs(reaction_rates1), axis=0)

    # Determine the maximum rate for each species-species pair and the maximum overall species-species rate
    max_species_rates1 = np.max(np.abs(species_rates1), axis=0)
    max_species_rate1 = np.max(max_species_rates1)
    species_index1 = max_species_rates1.reshape((num_species1 * num_species1)).argsort()
    max_species_rate_total = max_species_rate1
    max_concentration_total = max_concentration1


    nodes1 = []
    edges1 = []
    for i in range(num_species1 * num_species1):
        product_index1, reactant_index1 = divmod(species_index1[-i - 1], num_species1)
        if reactant_index1 > product_index1:
            # Both reactant -> product and product -> reactant are in this list, so only keep one of them
            continue
        if max_species_rates1[reactant_index1, product_index1] == 0:
            break
        if reactant_index1 not in nodes1 and len(nodes1) < max_node_count: nodes1.append(reactant_index1)
        if product_index1 not in nodes1 and len(nodes1) < max_node_count: nodes1.append(product_index1)
        if [reactant_index1, product_index1] not in edges1 and [product_index1, reactant_index1] not in edges1:
            edges1.append([reactant_index1, product_index1])
        if len(nodes1) > max_node_count:
            break
        if len(edges1) >= max_edge_count:
            break

    # Function to grab and generate the image for a given species
    def get_image_path(species):
        species_index = str(species) + '.png'
        image_path = ''
        if not species_path or not os.path.exists(species_path):  # species_path is defined while loading the mechanism
            raise OSError
        for root, dirs, files in os.walk(species_path):
            for f in files:
                if f == species_index:
                    image_path = os.path.join(root, f)
                    break
        if not image_path:
            image_path = os.path.join(species_path, species_index)
            species.molecule[0].draw(image_path)
        return image_path
    

    # Generate the flux diagram
    slope = -max_node_pen_width / np.log10(concentration_tol)

    # aramco_color = matplotlib.colors.to_hex((0.18627451, 0.48823529, 0.94117647))
    aramco_color = 'black'
    SAMPLE_RATE = 1
    for t1 in [a for a in range(len(distances1))][::SAMPLE_RATE]:  # time index
        graph = pydot.Dot('flux_diagram', graph_type='digraph', overlap="false")
        graph.set_fontname('sans')
        graph.set_fontsize('10')
        # ----------------------------ADD NODES ------------------------------#
        # For Mechanism 1
        for index1 in nodes1:
            nodewidths = np.zeros(3)  # keep track of species concentrations/nodewidths for all 3 mechanisms
            species1 = species_list1[index1]
            node1 = pydot.Node(name=str(species1))
            concentration1 = concentrations1[t1, index1] / max_concentration_total
            if concentration1 < concentration_tol:
                penwidth = 0.0
            else:
                penwidth = round(slope * np.log10(concentration1) + max_node_pen_width, 3)
                nodewidths[0] = penwidth
            node1.set_penwidth(penwidth)
            node1.set_fillcolor('white')
    #         node1.set_color(aramco_color)
            image_path1 = get_image_path(species1)
            if os.path.exists(image_path1):
                node1.set_image(image_path1)
                node1.set_label("")


            graph.add_node(node1)


        # ------------------------------- EDGES ------------------------------#
        # Add an edge for each species-species rate
        slope = -max_edge_pen_width / np.log10(species_rate_tol)

        # Go through edges in Mechanism 1
        for reactant_index1, product_index1 in edges1:
            if reactant_index1 in nodes1 and product_index1 in nodes1:
                reactant1 = species_list1[reactant_index1]
                product1 = species_list1[product_index1]
                label1 = ''
        #         label1 = label_flux_pair(reactant1, product1)

                edge1 = pydot.Edge(str(reactant1), str(product1), color=aramco_color)
                species_rate1 = species_rates1[t1, reactant_index1, product_index1] / max_species_rate_total
                if species_rate1 < 0:
                    edge1.set_dir("back")
                    species_rate1 = -species_rate1
                else:
                    edge1.set_dir("forward")
                # Set the edge pen width
                if species_rate1 < species_rate_tol:
                    penwidth = 0.0
                    edge1.set_dir("none")
                else:
                    penwidth = round(slope * np.log10(species_rate1) + max_edge_pen_width, 3)
                edge1.set_penwidth(penwidth)
                if label1 and label_automatically:
                    edge1.set_decorate(True)
                    edge1.set_label(label1)

                graph.add_edge(edge1)




        # General purpose graph settings
        graph.set_nodesep(0.11)
        graph.set_ranksep(0.35)
        graph.set_rankdir('LR')

        # Add Legend
        graph.add_node(pydot.Node(f't={distances1[t1]:.4e}', label=f't={distances1[t1]:.2e}s', color=aramco_color, shape='box', penwidth=max_node_pen_width))
        
    #     graph.add_node(pydot.Node(mech_1_label + f'\nt={times1[t1]:.4e}', label=mech_1_label + f'\nt={times1[t1]:.4e}', color=aramco_color, shape='box', penwidth=max_node_pen_width))
        # graph.add_node(pydot.Node(mech_2_label + f'\nt={times2[t2]:.4e}', label=mech_2_label + f'\nt={times2[t2]:.4e}', color=rmg1_color, shape='box', penwidth=max_node_pen_width))


        # write in multiple formats
    #     graph.write_dot(os.path.join(diagram_base_name, f'{diagram_base_name}_{t1:04}.dot')) # Yes this is supposed to be an index, not an actual time
        graph.write_png(os.path.join(working_dir, f'{diagram_base_name}_{t1:04}.png'))
    #     graph.write_pdf(os.path.join(diagram_base_name, f'{diagram_base_name}_{t1}.pdf'))
    #     graph.write_pdf(os.path.join(diagram_base_name, f'{diagram_base_name}_{t1}.svg'))
        if t1 > 40:  # don't generate too many
            break