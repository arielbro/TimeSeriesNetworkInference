import cProfile
import random
import time
import graphs
import attractors
from shutil import copyfile
# from pebble import ProcessPool
from collections import namedtuple
import csv
import stochastic
import multiprocessing
import itertools
import numpy
import os
import sys
import signal, os, errno
from functools import wraps
import stat
import platform

attractor_estimation_n_iter = 500
node_measurement_prob = 0.8
n_bio_experiments = 50
relax_experiments = False
timeout_seconds = int(0.3 * 60 * 60)
graph_parent_dir = "cellcollective_models"
max_attractor_length = 5
graph_size_filter = 10**6
tests_per_graph = 1
limit_experiments_to_short_attractors = False


def perform_graph_learning_tests(G, node_measurement_prob=node_measurement_prob, n_experiments=n_bio_experiments,
                                 relax_experiments=relax_experiments, timeout_seconds=timeout_seconds,
                                 max_attractor_length=max_attractor_length, tests_per_graph=tests_per_graph,
                                 limit_experiments_to_short_attractors=limit_experiments_to_short_attractors):
    """
    Given a model G, performs tests_per_graph tests, in each a random node's function is hidden, and a set
    of states in attractors are supplied to an ILP trying to find the function maximizing agreement with
    the states.
    :return: the average similarity between the learned function and hidden one in tests, the average proportion of
    experiment agreement in tests, and the average time taken per test.
    """
    # TODO: differentiate between tests in which the correct model was the correct solution vs otherwise.
    # TODO: better documentation here...

    test_times = []
    test_agreements = []
    test_model_similarities = []

    for _ in range(tests_per_graph):
        test_start = time.time() # TODO: measure ILP and experiment preparation time separately

        # hide a node's function
        G_hidden = G.copy()
        hidden_node_index = random.choice([i for i in range(len(G.vertices)) if
                                           len(G.vertices[i].predecessors()) > 0])
        G_hidden.vertices[hidden_node_index].function = None

        # create experiment measurements
        graph_attractors = stochastic.estimate_attractors(G, n_walks=attractor_estimation_n_iter, with_basins=False)
        if limit_experiments_to_short_attractors:
            graph_attractors = [att for att in graph_attractors if len(att) <= max_attractor_length]
            if len(graph_attractors) == 0:
                print "warning! No attractors for of suitable length found in graph. Experiments cannot reconstruct"
                return 0, 0, 0
        attractor_states = list(itertools.chain.from_iterable(graph_attractors))
        experiments = []
        for _ in range(n_experiments):
            experiment_state = random.choice(attractor_states)
            experiment = dict()
            for i in range(len(G.vertices)):
                # always measure input nodes
                if (len(G.vertices[i].predecessors()) == 0) or (random.random() < node_measurement_prob):
                    experiment[i] = experiment_state[i]
            assert len(experiment) != 0
            if len(experiment) != 0:
                experiments.append(experiment)
            else:
                print "Warning: Experiment came out with no measured nodes"


        # try to learn the model back
        G_found, agreement = attractors.learn_model_from_experiment_agreement(
            G_hidden, experiments, relax_experiments,
            max_attractor_length, timeout_seconds)
        model_similarity = 0
        for var_combination in itertools.product((False, True),
                                                 repeat=len(G.vertices[hidden_node_index].predecessors())):
            if G.vertices[hidden_node_index].function(*var_combination) == \
                    G_found.vertices[hidden_node_index].function(*var_combination):
                model_similarity += 1
        model_similarity /= float(2 ** len(G.vertices[hidden_node_index].predecessors()))

        test_model_similarities.append(model_similarity)
        test_agreements.append(agreement)
        test_times.append(time.time() - test_start)
    average_similarity = sum(test_model_similarities) / float(len(test_model_similarities))
    average_test_agreement = sum(test_agreements) / float(len(test_agreements))
    average_test_time = sum(test_times) / float(len(test_times))
    return average_similarity, average_test_agreement, average_test_time


ModelLearningResult = namedtuple("ModelLearningResult", "graph_name size "
                                                    "n_inputs normalized_n_inputs "
                                                    "max_degree mean_degree "
                                                    "node_measurement_prob n_bio_experiments "
                                                    "relax_experiments max_attractor_length "
                                                    "limit_experiments_to_short_attractors "
                                                    "average_model_similarity "
                                                    "average_experiment_agreement "
                                                    "average_test_time"
                                )


def main():
    print ("node_measurement_prob={}, n_experiments={}, relax_experiments={}, timeout_seconds={}, "
           "graph_parent_dir={}, max_attractor_length={}, graph_size_filter={}, n_tests={}".format(
            node_measurement_prob, n_bio_experiments, relax_experiments, timeout_seconds, graph_parent_dir,
            max_attractor_length, graph_size_filter, tests_per_graph))

    # TODO: use grownups' argument parsing library.
    if len(sys.argv) == 1:
        output_path = 'model_learning_statistics.csv'
    else:
        output_path = sys.argv[1]
    print "saving output to {}".format(output_path)

    # filter graphs
    tested_graphs = []
    graph_names = []
    candidate_graph_names = os.listdir(graph_parent_dir)
    filter_start = time.time()
    for graph_dir in candidate_graph_names:
        try:
            G = graphs.Network.parse_boolean_tables(os.path.join(graph_parent_dir, graph_dir))
            if len(G.vertices) <= graph_size_filter:
                tested_graphs.append(G)
                graph_names.append(graph_dir)

        except ValueError as e:
            if e.message.startswith("Model export from cellcollective failed"):
                print "warning - did not load graph {}".format(graph_dir)

    print "Filtering done, time taken: {:.2f} secs. {} graphs left".format(time.time() - filter_start,
                                                                           len(tested_graphs))

    # perform tests
    results = []
    for graph, graph_name in zip(tested_graphs, graph_names):
        graph_start = time.time()
        print "testing graph {}".format(graph_name)

        average_similarity, average_test_agreement, average_test_time = perform_graph_learning_tests(graph)

        # get general attributes of graph
        n_inputs = len([v for v in graph.vertices if len(v.predecessors()) == 0])
        max_degree = max([len(v.predecessors()) for v in graph.vertices])
        size = len(graph.vertices)
        mean_degree = sum([len(v.predecessors()) for v in graph.vertices]) / float(size)
        normaliezd_n_inputs = n_inputs / float(size)

        result = ModelLearningResult(graph_name, size, n_inputs, normaliezd_n_inputs, max_degree, mean_degree,
                                     node_measurement_prob, n_bio_experiments, relax_experiments, max_attractor_length,
                                     average_similarity, average_test_agreement, average_test_time)
        results.append(result)
        print "time taken for graph {}: {:.2f} secs".format(graph_name, time.time() - graph_start)

        # Re-save results after each graph is tested.
        with open(output_path, 'wb') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["graph_name", "size", "n_inputs", "normalized_n_inputs", "max_degree", "mean_degree",
                             "node_measurement_prob", "n_bio_experiments", "relax_experiments", "max_attractor_length",
                             "average_model_similarity", "average_experiment_agreement", "average_test_time",
                             "limit_experiments_to_short_attractors"
                             ])
            for result in results:
                if result is None:
                    continue
                writer.writerow(result)


if __name__ == "__main__":
    main()
