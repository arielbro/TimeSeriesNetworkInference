import os
import shutil
from attractor_learning import graphs
import numpy as np
import itertools
import logging
import random


def generate_random_graphs(graphs_dir, **kwargs):
    reference_graphs = []
    for graph_dir in os.listdir(graphs_dir):
        G = graphs.Network.parse_boolean_tables(os.path.join(graphs_dir, graph_dir))
        reference_graphs.append(G)

    for G in reference_graphs:
        for _ in range(kwargs['random_networks_per_reference']):
            random_graph = G.copy()
            random_graph.randomize_edges_by_switching(n_attempts=1000)
            random_graph.randomize_functions(mutate_input_nodes=kwargs['mutate_input_nodes'],
                                             preserve_truth_ratio=kwargs['preserve_truth_ratio'],
                                             function_type_restriction=kwargs['function_type_restriction'])
            assert random_graph != G
            yield random_graph


def generate_scaffold_network(G, **kwargs):
    """
    Generates an almost correct scaffold network to use in model inference.
    Currently just adds new edges to the reference graph G.
    :param added_edge_frac: fraction of current amount of edges to add. |E'|=(1+added_edge_frac)|E|
    :param G:
    :return:
    """
    added_edge_frac = kwargs['scaffold_network_added_edge_fraction']
    removed_edge_frac = kwargs['scaffold_network_removed_edge_fraction']
    preserve_input_nodes_on_add = kwargs['preserve_input_nodes_on_add']

    scaffold = G.copy()
    for vertex in scaffold.vertices:
        vertex.function = None

    n_removed_edges = int(len(scaffold.edges) * removed_edge_frac)
    removed_edges = random.sample(scaffold.edges, n_removed_edges)

    n_added_edges = int(len(scaffold.edges) * added_edge_frac)
    optional_added_edges = list((a, b) for (a, b) in itertools.combinations(scaffold.vertices, 2) if
                      (a, b) not in scaffold.edges)
    if preserve_input_nodes_on_add:
        optional_added_edges = [(a, b) for (a, b) in optional_added_edges if len(b.predecessors()) > 0]
    if n_added_edges > len(optional_added_edges):
        warning = "Warning! More edges to add than possible. Reducing amount of added edges"
        logger = logging.getLogger(__name__)
        logger.warning(warning)
        print(warning)
        n_added_edges = len(optional_added_edges)
    added_edges = random.sample(optional_added_edges, n_added_edges)  # without replacement
    # TODO: edge (haha) cases (e.g. not enough edges to choose from)

    for edge in removed_edges:
        scaffold.edges.remove(edge)
    scaffold.edges.extend(added_edges)
    for node in scaffold.vertices:
        node.precomputed_predecessors = None
        node.precomputed_successors = None
    return scaffold
