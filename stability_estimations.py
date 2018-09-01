import random
import time
import graphs
import attractors
from collections import namedtuple
import csv
import numpy
import os
import sys

if __name__ == "__main__":
    # TODO: use grownups' argument parsing library.
    if len(sys.argv) == 1:
        output_path = 'bitchange_results.csv'
        analyze_originals = True
    else:
        output_path = sys.argv[1]
        analyze_originals = bool(int(sys.argv[2]))

    experiment_max_len = 15
    experiment_n_iter = 300

    StabilityResult = namedtuple("StabilityResult", "graph_name random_functions random_edges size "
                                                    "n_inputs normalized_n_inputs max_degree mean_degree"
                                                    "state_bits_changed"
                                                    # "minimal_model_bitchange "
                                                    "model_bitchange_prob "
                                                    "state_bitchange_prob"
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
        graph_name_to_attributes[graph.name] = dict(n_inputs=n_inputs, max_degree=max_degree,
                                                    size=size, mean_degree=mean_degree,
                                                    normaliezd_n_inputs=normaliezd_n_inputs)
        print "#{}; {} input nodes for graph {} of size {} and max degree {}".format(i, n_inputs, name,
                                                                                     size, max_degree)
    if analyze_originals:
        for i, graph, name in zip(range(len(biological_graphs)), biological_graphs, biological_graph_names):
            print "working on biological graph #{}".format(i)
            start = time.time()
            # try:
            #     minimal_model_bitchange = \
            #         attractors.find_model_bitchange_for_new_attractor(graph, max_len=experiment_max_len,
            #                                                           verbose=False, use_dubrova=False,
            #                                                           simplify_boolean=False)
            # except attractors.TimeoutError as e:
            #     minimal_model_bitchange = numpy.inf
            print "time taken for model_bitchange={:.2f} secs".format(time.time() - start)
            second_start = time.time()
            state_bits_changed = random.randint(5)
            try:
                model_bitchange_probability = \
                    attractors.find_model_bitchange_probability_for_different_attractors(graph,
                                                                                         n_iter=experiment_n_iter,
                                                                                         use_dubrova=False)
            except attractors.TimeoutError as e:
                model_bitchange_probability = numpy.inf
            print "time taken for model_bitchange_probability={:.2f} secs".format(time.time() - second_start)
            third_start = time.time()
            try:
                state_bitchange_probability = \
                    attractors.find_state_bitchange_probability_for_different_attractors(graph,
                                                                                         n_iter=experiment_n_iter,
                                                                                         n_bits=state_bits_changed,
                                                                                         parallel=True)
            except attractors.TimeoutError as e:
                state_bitchange_probability = numpy.inf
            print "time taken for state_bitchange_probability={:.2f} secs".format(time.time() - third_start)

            result = StabilityResult(graph_name=name, random_functions=False, random_edges=False,
                                     size=graph_name_to_attributes[name]['size'],
                                     n_inputs=graph_name_to_attributes[name]['n_inputs'],
                                     normalized_n_inputs=graph_name_to_attributes[name]['normalized_n_inputs'],
                                     max_degree=graph_name_to_attributes[name]['max_degree'],
                                     mean_degree=graph_name_to_attributes[name]['mean_degree'],
                                     # minimal_model_bitchange=minimal_model_bitchange,
                                     model_bitchange_prob=model_bitchange_probability,
                                     state_bits_changed=state_bits_changed,
                                     state_bitchange_prob=state_bitchange_probability
                                     )
            results.append(result)

            # print "bitchange={}".format(minimal_model_bitchange)
            # print "model_prob={:.2f}, state_prob={:.2f}".format(
                # minimal_model_bitchange,
                # model_bitchange_probability,
                # state_bitchange_probability)
            print "time_taken={:.2f} secs".format(time.time() - start)

    for test in range(5000):
        print "test #{}".format(test + 1)
        for graph, name in zip(biological_graphs, biological_graph_names):
            graph_copy = graph.copy()

            randomize_functions, randomize_edges = random.choice([(False, True), (True, False), (True, True)])
            if randomize_functions:
                graph_copy.randomize_functions(preserve_truth_ratio=True)
            if randomize_edges:
                graph_copy.randomize_edges()
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
            try:
                model_bitchange_probability = \
                    attractors.find_model_bitchange_probability_for_different_attractors(graph_copy,
                                                                                         max_len=experiment_max_len,
                                                                                         n_iter=experiment_n_iter,
                                                                                         use_dubrova=False)
            except attractors.TimeoutError as e:
                model_bitchange_probability = numpy.inf
            print "time taken for model_bitchange_probability={:.2f} secs".format(time.time() - second_start)
            third_start = time.time()
            try:
                state_bitchange_probability = \
                    attractors.find_state_bitchange_probability_for_different_attractors(graph_copy,
                                                                                         n_iter=experiment_n_iter,
                                                                                         n_bits=state_bits_changed,
                                                                                         parallel=True)
            except attractors.TimeoutError as e:
                state_bitchange_probability = numpy.inf
            print "time taken for state_bitchange_probability={:.2f} secs".format(time.time() - third_start)

            result = StabilityResult(graph_name=name, random_functions=randomize_functions,
                                     random_edges=randomize_edges,
                                     size=graph_name_to_attributes[name]['size'],
                                     n_inputs=graph_name_to_attributes[name]['n_inputs'],
                                     normalized_n_inputs=graph_name_to_attributes[name]['normalized_n_inputs'],
                                     max_degree=graph_name_to_attributes[name]['max_degree'],
                                     mean_degree=graph_name_to_attributes[name]['mean_degree'],
                                     # minimal_model_bitchange=minimal_model_bitchange,
                                     model_bitchange_prob=model_bitchange_probability,
                                     state_bits_changed=state_bits_changed,
                                     state_bitchange_prob=state_bitchange_probability
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
                             "minimal_model_bitchange",
                             # "model_bitchange_prob",
                             # "state_bitchange_prob"
                             ])
            for stability_result in results:
                writer.writerow(stability_result)