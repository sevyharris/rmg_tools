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

import cantera as ct
import pydot


def get_i_thing(ref_composition, phase):
    for i in range(phase.n_species):
        if phase.species()[i].composition == ref_composition:
            return i
    assert False


gas = ct.Solution('/home/moon/rmg/RMG-Py/examples/rmg/ethane-oxidation/cantera/chem_annotated.yaml')


i_ethane = get_i_thing({'C': 2, 'H':6}, gas)

gas.species()


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
graph = pydot.Dot('flux_diagram', graph_type='digraph', overlap="false")

for node in nodes:
    pydot_node = pydot.Node(name=node.name)
    # node1.set_penwidth(penwidth)
    # node1.set_fillcolor('white')
    # node1.set_color(aramco_color)
    # node1.set_image(image_path1)
    # node1.set_label("")
    graph.add_node(pydot_node)

for edge in edges:
    pydot_edge = pydot.Edge(edge.r_node.name, edge.p_node.name)

    # edge1.set_dir("back")
    # edge1.set_dir("forward")
    # edge1.set_dir("none")
    # edge1.set_decorate(True)
    # edge1.set_label(label1)

    graph.add_edge(pydot_edge)

# # General purpose graph settings
# graph.set_nodesep(0.11)
# graph.set_ranksep(0.35)
graph.set_rankdir('LR')
graph.write_png('example.png')
# -



gas.reactions()[0].reactants





start_species = gas.species_names[i_ethane]

for 



start_species


