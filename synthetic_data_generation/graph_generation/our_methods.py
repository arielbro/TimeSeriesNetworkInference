import os
import shutil
from attractor_learning import graphs
import numpy as np
import itertools
import logging
import random


# scaffold_network_added_edge_fraction value asking for the full scaffold rather than a fraction of added
# edges. Kept here rather than in generate.py so the option parser and this function agree on the spelling.
FULL_SCAFFOLD = "all"


def parse_added_edge_fraction(value):
    """configargparse `type` for --scaffold_network_added_edge_fraction: a float, or FULL_SCAFFOLD ("all")
    for a scaffold holding every possible edge (see generate_scaffold_network)."""
    if str(value).strip().lower() == FULL_SCAFFOLD:
        return FULL_SCAFFOLD
    return float(value)


def generate_random_graphs(graphs_dir, **kwargs):
    reference_graphs = []
    for graph_dir in os.listdir(graphs_dir):
        # a graphs_dir may hold loose files next to the model directories (e.g. a README)
        if not os.path.isdir(os.path.join(graphs_dir, graph_dir)):
            continue
        # max_graph_size bounds the number of nodes; None (the default, i.e. the option left out of the
        # config) means no filtering. Probed from the model's index/header so an oversized model is never
        # parsed, in whichever of the two model-directory formats it is stored.
        if (kwargs['max_graph_size'] is not None) and \
                (graphs.Network.model_dir_size(os.path.join(graphs_dir, graph_dir))
                 > kwargs['max_graph_size']):
            continue
        print("loading graph {}".format(graph_dir))
        G = graphs.Network.parse_model_dir(os.path.join(graphs_dir, graph_dir))
        reference_graphs.append(G)

    for G in reference_graphs:
        for _ in range(kwargs['random_networks_per_reference']):
            yield randomize_reference_graph(G, **kwargs)


def randomize_reference_graph(reference, **kwargs):
    """One randomization of a reference graph: degree-preserving edge switching, then fresh functions.
    Split out of generate_random_graphs so a worker process can randomize a reference it parsed itself -
    a Network holding truth-table functions cannot be pickled, so it cannot be handed to one."""
    random_graph = reference.copy()
    random_graph.randomize_edges_by_switching(n_attempts=1000)
    random_graph.randomize_functions(mutate_input_nodes=kwargs['mutate_input_nodes'],
                                     preserve_truth_ratio=kwargs['preserve_truth_ratio'],
                                     function_type_restriction=kwargs['function_type_restriction'])
    assert random_graph != reference
    return random_graph


def generate_scaffold_network(G, **kwargs):
    """
    Generates an almost correct scaffold network to use in model inference.
    Currently just adds new edges to the reference graph G.
    :param added_edge_frac: fraction of current amount of edges to add. |E'|=(1+added_edge_frac)|E|.
        The string "all" instead gives the full scaffold - every possible edge, so inference is
        unconstrained by the scaffold's topology - in which case removed_edge_frac is ignored.
    :param G:
    :return:
    """
    added_edge_frac = kwargs['scaffold_network_added_edge_fraction']
    removed_edge_frac = kwargs['scaffold_network_removed_edge_fraction']
    preserve_input_nodes_on_add = kwargs['preserve_input_nodes_on_add']

    scaffold = G.copy()
    for vertex in scaffold.vertices:
        vertex.function = None

    if added_edge_frac == FULL_SCAFFOLD:
        # every ordered pair, self-loops included, so the scaffold is a superset of G's edges whatever they
        # are. Removing edges would contradict that, so removed_edge_frac has no meaning here and is ignored.
        # Targets are picked before the edges are replaced, since predecessors() then still describes G.
        targets = [v for v in scaffold.vertices
                   if (not preserve_input_nodes_on_add) or len(v.predecessors()) > 0]
        scaffold.edges = [(u, v) for v in targets for u in scaffold.vertices]
        for node in scaffold.vertices:
            node.precomputed_predecessors = None
            node.precomputed_successors = None
        return scaffold

    n_removed_edges = int(len(scaffold.edges) * removed_edge_frac)
    removed_edges = random.sample(scaffold.edges, n_removed_edges)

    n_added_edges = int(len(scaffold.edges) * added_edge_frac)
    # membership against a set: scaffold.edges is a list, and an `in` test per candidate pair is
    # O(#edges), which on a 1000-node graph turns this into ~10**9 comparisons
    existing_edges = set(scaffold.edges)
    # every ordered pair the graph doesn't already have, self-loops included - the same candidate set the
    # FULL_SCAFFOLD branch lays down, just sampled from. itertools.combinations, used here previously, yields
    # each pair once in vertex order, so only about half the missing edges could ever be drawn and every
    # added edge pointed from a lower-indexed node to a higher-indexed one.
    optional_added_edges = [(a, b) for (a, b) in itertools.product(scaffold.vertices, repeat=2)
                            if (a, b) not in existing_edges]
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
