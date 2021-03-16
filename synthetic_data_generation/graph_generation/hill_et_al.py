import networkx
import numpy as np
import random
import itertools

"""
Implementing (if I got them right) the simulations from
https://academic.oup.com/bioinformatics/article/28/21/2804/235527
"""


def generate_random_graph(n_nodes=20, n_edges=74):
    """
    Generate a random graph based on an Erdos-Renyi model.
    :param n_nodes:
    :param n_edges:
    :return:
    """
    return networkx.generators.gnm_random_graph(n=n_nodes, m=n_edges, directed=True)


def prior_graph_from_data_graph(data_graph, n_removed=54, n_added=10):
    """
    Create a random graph from a given one, to be used as a prior graph for inference.
    Done by selecting n_removed edges at uniform and removing them, and n_added edges
    not in the original graph to add, at uniform.
    :param data_graph:
    :param n_removed:
    :param n_added:
    :return:
    """
    G = data_graph.copy()
    assert n_removed <= len(G.edges)
    assert n_added <= len(G.nodes) * (len(G.nodes) - 1) - len(G.edges)
    removed_edges = random.sample(G.edges, n_removed)
    added_edges = random.sample([e for e in itertools.permutations(G.edges, 2) if
                                 e not in G.edges], n_added)
    for e in removed_edges:
        G.remove_edge(*e)
    for e in added_edges:
        G.add_edge(*e)
    return G

