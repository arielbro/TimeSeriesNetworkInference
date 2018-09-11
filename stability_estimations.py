import random
import time
import graphs
import attractors
from shutil import copyfile
from collections import namedtuple
import csv
import numpy
from stochastic import estimate_probability_to_reroll_attractor
import os
import sys
import stat

if __name__ == "__main__":
    # TODO: use grownups' argument parsing library.
    if len(sys.argv) == 1:
        output_path = 'bitchange_results.csv'
    else:
        output_path = sys.argv[1]

    experiment_max_len = 15
    experiment_n_iter = 100
    biological_test_frequency = 10

    # copy Dubrova's executable, to prevent mysterious errors when running multiple processes calling it.
    process_specific_dubrova_path = "temp_{}_{}".format(attractors.dubrova_path, os.getpid())
    copyfile(attractors.dubrova_path, process_specific_dubrova_path)
    os.chmod(process_specific_dubrova_path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
    attractors.dubrova_path = process_specific_dubrova_path

    StabilityResult = namedtuple("StabilityResult", "graph_name random_functions random_edges size "
                                                    "n_inputs normalized_n_inputs max_degree mean_degree "
                                                    "state_bits_changed "
                                                    # "minimal_model_bitchange "
                                                    "model_bitchange_prob "
                                                    "state_bitchange_prob "
                                                    "reroll_prob"
                                )
    results = []

    biological_graphs = []
    biological_graph_names = os.listdir("cellcollective_models")
    for graph_dir in biological_graph_names:
        try:
            G = graphs.Network.parse_boolean_tables(os.path.join("cellcollective_models", graph_dir))
            biological_graphs.append(G)
        except ValueError as e:
            if e.message.startswith("Model export from cellcollective failed"):
                print "warning - did not load graph {}".format(graph_dir)

    graph_name_to_attributes = dict()
    for i, graph, name in zip(range(len(biological_graphs)), biological_graphs, biological_graph_names):
        n_inputs = len([v for v in graph.vertices if len(v.predecessors()) == 0])
        max_degree = max([len(v.predecessors()) for v in graph.vertices])
        size = len(graph.vertices)
        mean_degree = sum([len(v.predecessors()) for v in graph.vertices]) / float(size)
        normaliezd_n_inputs = n_inputs / float(size)
        graph_name_to_attributes[name] = dict(n_inputs=n_inputs, max_degree=max_degree,
                                                    size=size, mean_degree=mean_degree,
                                              normalized_n_inputs=normaliezd_n_inputs)
        print "#{}; {} input nodes for graph {} of size {} and max degree {}".format(i, n_inputs, name,
                                                                                     size, max_degree)
    for test in range(5000):
        print "test #{}".format(test)
        for graph, name in zip(biological_graphs, biological_graph_names):
            is_biological = (test % biological_test_frequency) == 0
            graph_copy = graph.copy()

            if not is_biological:
                randomize_functions, randomize_edges = random.choice([(False, True), (True, False), (True, True)])
                if randomize_functions:
                    graph_copy.randomize_functions(preserve_truth_ratio=True)
                if randomize_edges:
                    # graph_copy.randomize_incoming_edges()
                    graph_copy.randomize_edges_by_switching()
            start = time.time()

            # try:
            #     minimal_model_bitchange = \
            #         attractors.find_model_bitchange_for_new_attractor(graph_copy, max_len=experiment_max_len,
            #                                                           verbose=False, use_dubrova=False,
            #                                                           simplify_boolean=False)
            # except attractors.TimeoutError as e:
            #     minimal_model_bitchange = numpy.inf
            print "time taken for model_bitchange={:.2f} secs".format(time.time() - start)
            second_start = time.time()
            state_bits_changed = random.randint(1, 5)
            try:
                model_bitchange_probability = \
                    attractors.find_model_bitchange_probability_for_different_attractors(graph_copy,
                                                                                         max_len=experiment_max_len,
                                                                                         n_iter=experiment_n_iter,
                                                                                         use_dubrova=True)
            except attractors.TimeoutError as e:
                model_bitchange_probability = numpy.inf
            print "time taken for model_bitchange_probability={:.2f} secs".format(time.time() - second_start)
            third_start = time.time()
            try:
                state_bitchange_probability = \
                    attractors.find_state_bitchange_probability_for_different_attractors(graph_copy,
                                                                                         n_iter=experiment_n_iter,
                                                                                         n_bits=state_bits_changed,
                                                                                         parallel=False)
            except attractors.TimeoutError as e:
                state_bitchange_probability = numpy.inf
            print "time taken for state_bitchange_probability={:.2f} secs".format(time.time() - third_start)

            fourth_start = time.time()
            reroll_probability = estimate_probability_to_reroll_attractor(graph_copy, n_walks=experiment_n_iter)
            print "time taken for reroll probability={:.2f} secs".format(time.time() - fourth_start)

            result = StabilityResult(graph_name=name,
                                     random_functions=randomize_functions if not is_biological else False,
                                     random_edges=randomize_edges if not is_biological else False,
                                     size=graph_name_to_attributes[name]['size'],
                                     n_inputs=graph_name_to_attributes[name]['n_inputs'],
                                     normalized_n_inputs=graph_name_to_attributes[name]['normalized_n_inputs'],
                                     max_degree=graph_name_to_attributes[name]['max_degree'],
                                     mean_degree=graph_name_to_attributes[name]['mean_degree'],
                                     # minimal_model_bitchange=minimal_model_bitchange,
                                     model_bitchange_prob=model_bitchange_probability,
                                     state_bits_changed=state_bits_changed,
                                     state_bitchange_prob=state_bitchange_probability,
                                     reroll_prob=reroll_probability
                                     )

            # print "model bitchange={}".format(minimal_model_bitchange)
            results.append(result)
            print "model_prob={:.2f}, state_prob={:.2f}".format(
                # minimal_model_bitchange,
                model_bitchange_probability,
                state_bitchange_probability
            )
            print "time_taken = {:.2f} secs".format(time.time() - start)

        # save on each iteration, why not
        with open(output_path, 'wb') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["graph_name", "random_functions", "random_edges", "n",
                             "n_inputs", "normalized_n_inputs", "max_degree", "mean_degree",
                             "state_bits_changed",
                             "model_bitchange_prob",
                             "state_bitchange_prob",
                             "reroll_prob"
                             ])
            for stability_result in results:
                writer.writerow(stability_result)
