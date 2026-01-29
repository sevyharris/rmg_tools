# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import cantera as ct
import pydot
import os
import numpy as np


# -

def get_i_thing(ref_composition, phase):
    for i in range(phase.n_species):
        if phase.species()[i].composition == ref_composition:
            return i
    assert False


gas = ct.Solution('/home/moon/rmg/RMG-Py/examples/rmg/ethane-oxidation/cantera/chem_annotated.yaml')


i_ethane = get_i_thing({'C': 2, 'H':6}, gas)
i_Ar = get_i_thing({'Ar': 1}, gas)
i_O2 = get_i_thing({'O': 2}, gas)

# +
# run a shock tube simulation to get flux

T = 800  # K
P = ct.one_atm
X = {gas.species_names[i_Ar]: 0.7649, gas.species_names[i_ethane]: 0.03135, gas.species_names[i_O2]: 0.2038}
t_end = 1.0  # s

gas.TPX = T, P, X
reactor = ct.IdealGasReactor(gas)
reactor_net = ct.ReactorNet([reactor])

times = [0]
T = [reactor.T]
P = [reactor.thermo.P]
X = [reactor.thermo.X]  # mol fractions
rates = [reactor.kinetics.net_rates_of_progress]
while reactor_net.time < t_end:
    reactor_net.step()
    
    times.append(reactor_net.time)
    T.append(reactor.T)
    P.append(reactor.thermo.P)
    X.append(reactor.thermo.X)
    rates.append(reactor.kinetics.net_rates_of_progress)

times = np.array(times)
P = np.array(P)
X = np.array(X)
rates = np.array(rates)
# -

X.shape

rates.shape

# +
# at any given time, give me 
# -



# +
# build a graph of the reaction network

class Node():
    def __init__(self, name, index):
        self.name = name
        self.index = index
        self.edges = []
    def add_edge(self, edge):
        edge_indices = [e.index for e in self.edges]
        if edge.index not in edge_indices:
            self.edges.append(edge)


class Edge():
    # edges with the same name and index can connect different sets of nodes
    # this is my solution to not having flux pairs
    def __init__(self, name, index, r_node, p_node):
        self.name = name
        self.index = index
        self.r_node = r_node
        self.p_node = p_node
        self.r_node.add_edge(self)
        self.p_node.add_edge(self)
    def __repr__(self):
        return self.name


# +
nodes = []
for i in range(gas.n_species):
    nodes.append(Node(gas.species_names[i], i))

def get_node_index_by_name(name):
    return gas.species_names.index(name)

edges = []
for i in range(gas.n_reactions):
    for r in gas.reactions()[i].reactants.keys():
        for p in gas.reactions()[i].products.keys():
            # connect the edges together
            r_node = nodes[get_node_index_by_name(r)]
            p_node = nodes[get_node_index_by_name(p)]
            edges.append(Edge(gas.reactions()[i].equation, i, r_node, p_node))


# +
species_image_dir = '/home/moon/rmg/RMG-Py/examples/rmg/ethane-oxidation/chemkin/species_images/'
# make images bigger

# for i in range(gas.n_species):

#     cairosvg.svg2png(
#         url=os.path.join(species_image_dir, f'{gas.species_names[i]}.svg'),
#         write_to=os.path.join(species_image_dir, f'{gas.species_names[i]}.png'),
#         scale=1.2,
        
#     )


# -

SPEC_CONC_TOL

# +
graph = pydot.Dot('flux_diagram', graph_type='digraph', overlap="false")

RATE_FLUX_TOL = 1e-9
SPEC_CONC_TOL = 1e-8
MAX_SPEC_CONC = np.max(X)
MAX_FLUX = np.max(np.abs(rates))
MAX_NODE_PEN_WIDTH = 9
MAX_EDGE_PEN_WIDTH = 9
SLOPE_NODE = -MAX_NODE_PEN_WIDTH / np.log10(SPEC_CONC_TOL)
SLOPE_EDGE = -MAX_EDGE_PEN_WIDTH / np.log10(RATE_FLUX_TOL)

t = 100  # time snapshot

for t in [0, 20, 50, 100, 200]:
    
    for node in nodes:
        if node.edges == []:
            continue
    
        spec_conc = X[t, node.index] / MAX_SPEC_CONC
        if spec_conc < SPEC_CONC_TOL:
            node_penwidth = 0.0
        else:
            node_penwidth = round(SLOPE_NODE * np.log10(spec_conc) + MAX_NODE_PEN_WIDTH, 3)
    
        pydot_node = pydot.Node(name=node.name)  # , width=2.5, height=4
        
        # pydot_node.set_width(2.0)
        pydot_node.set_image(os.path.join(species_image_dir, f'{node.name}.png'))
        pydot_node.set_label('')
        pydot_node.set_fillcolor('white')
        pydot_node.set_color('black')
        pydot_node.set_penwidth(node_penwidth)
        
        graph.add_node(pydot_node)
    
    unique_edges = []
    edge_indices_plotted = []
    for edge in edges:
        # if edge.index in edge_indices_plotted:
        #     continue
        # edge_indices_plotted.append(edge.index)
        # unique_edges.append(edge)
        pydot_edge = pydot.Edge(edge.r_node.name, edge.p_node.name)
    
    
        flux = rates[t, edge.index]
        pydot_edge.set_dir("forward")
        if flux < 0:
            pydot_edge.set_dir("back")
    
        flux = np.abs(flux)
        if flux < RATE_FLUX_TOL:
            continue
            edge_penwidth = 0.0
            pydot_edge.set_dir("none")
        else:
            edge_penwidth = round(SLOPE_EDGE * np.log10(flux) + MAX_EDGE_PEN_WIDTH, 3)
    
        pydot_edge.set_penwidth(edge_penwidth)
        # edge1.set_dir("none")
        # edge1.set_decorate(True)
        # edge1.set_label(label1)
    
        graph.add_edge(pydot_edge)
    
    # # General purpose graph settings
    # graph.set_nodesep(0.11)
    # graph.set_ranksep(0.35)
    graph.set_rankdir('LR')
    graph.write_png(f'example_{t:04}.png')
    
    # assert len(unique_edges) == gas.n_reactions
# -



gas.reactions()[0].reactants





start_species = gas.species_names[i_ethane]

for 



start_species


