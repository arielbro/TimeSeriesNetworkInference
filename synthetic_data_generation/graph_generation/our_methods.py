import os
import shutil
from attractor_learning import graphs
import numpy as np
import itertools
import logging
import random
from attractor_learning.ilp import FunctionTypeRestriction

random_networks_per_reference = 10
graphs_dir = "../data/cellcollective_sample"
mutate_input_nodes = True
preserve_truth_ratio = True
function_type_restriction = FunctionTypeRestriction.NONE

scaffold_network_added_edge_fraction = 0.2

reference_graphs = []
for graph_dir in os.listdir(graphs_dir):
    G = graphs.Network.parse_boolean_tables(os.path.join(graphs_dir, graph_dir))
    reference_graphs.append(G)


def generate_random_graphs(logger=None):
    for G in reference_graphs:
        for _ in range(random_networks_per_reference):
            random_graph = G.copy()
            random_graph.randomize_edges_by_switching(n_attempts=1000)
            random_graph.randomize_functions(mutate_input_nodes=mutate_input_nodes,
                                             preserve_truth_ratio=preserve_truth_ratio,
                                             function_type_restriction=function_type_restriction)
            assert random_graph != G
            yield random_graph


def generate_scaffold_network(G, added_edge_frac=scaffold_network_added_edge_fraction, logger=None):
    """
    Generates an almost correct scaffold network to use in model inference.
    Currently just adds new edges to the reference graph G.
    :param added_edge_frac: fraction of current amount of edges to add. |E'|=(1+added_edge_frac)|E|
    :param G:
    :return:
    """
    scaffold = G.copy()
    for vertex in scaffold.vertices:
        vertex.function = None
    n_added_edges = int(len(scaffold.edges) * added_edge_frac)
    optional_edges = list((a, b) for (a, b) in itertools.combinations(scaffold.vertices, 2) if
                      (a, b) not in scaffold.edges)
    if n_added_edges > len(optional_edges):
        warning = "Warning! More edges to add than possible. Reducing amount of added edges"
        if logger is not None:
            logger.warning(warning)
        print(warning)
        n_added_edges = len(optional_edges)
    added_edges = random.sample(optional_edges, n_added_edges)  # without replacement
    # TODO: edge (haha) cases (e.g. not enough edges to choose from)
    scaffold.edges.extend(added_edges)
    for node in scaffold.vertices:
        node.precomputed_predecessors = None
        node.precomputed_successors = None
    return scaffold


def log_params():
    logger = logging.getLogger(__name__)
    logger.info("random_networks_per_reference={}".format(random_networks_per_reference))
    logger.info("graphs_dir={}".format(graphs_dir))
    logger.info("mutate_input_nodes={}".format(mutate_input_nodes))
    logger.info("preserve_truth_ratio={}".format(preserve_truth_ratio))
    logger.info("scaffold_network_added_edge_fraction={}".format(scaffold_network_added_edge_fraction))
    logger.info("function_type_restriction={}".format(function_type_restriction))