import random
import time
import graphs
import attractors
from collections import namedtuple
import csv
import numpy

experiment_max_len = 15
experiment_n_iter = 100

StabilityResult = namedtuple("StabilityResult", "graph_name random_functions random_edges n "
                                                "minimal_model_bitchange "
                                                "model_bitchange_prob "
                                                "state_bitchange_prob")

results = []
biological_graphs = [
    graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
                              "BNS_Dubrova_2011\\MAPK_large2.cnet"),
    # graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
    #                           "BNS_Dubrova_2011\\arabidopsis.cnet.txt"),
    # takes a longest to load (~14 secs), it stalls the experiments as well (has max degree of 14 !!!)
    # graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
    #                           "BNS_Dubrova_2011\\EGFR_man.cnet"),
    graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
                              "BNS_Dubrova_2011\\MAPK_large.cnet"),
    graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
                              "BNS_Dubrova_2011\\thelper.cnet.txt"),
    graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
                              "BNS_Dubrova_2011\\tcr.cnet"),
]
biological_graph_names = ["MAPK_large2", "arabidopsis", "EGFR_man", "MAPK_large", "thelper", "tcr"]

for graph, name in zip(biological_graphs, biological_graph_names):
    n_inputs = len([v for v in graph.vertices if len(v.predecessors()) == 0])
    max_degree = max([len(v.predecessors()) for v in graph.vertices])
    print "{} input nodes for graph {} of size {} and max degree {}".format(n_inputs, name, len(graph.vertices),
                                                                            max_degree)

for graph, name in zip(biological_graphs, biological_graph_names):
    start = time.time()
    try:
        minimal_model_bitchange = \
            attractors.find_model_bitchange_for_new_attractor(graph, max_len=experiment_max_len,
                                                              verbose=False, use_dubrova=True,
                                                              simplify_boolean=False)
    except attractors.TimeoutError as e:
        minimal_model_bitchange = numpy.inf
    print "time taken for model_bitchange={:.2f} secs".format(time.time() - start)
    second_start = time.time()
    try:
        model_bitchange_probability = \
            attractors.find_model_bitchange_probability_for_different_attractors(graph, n_iter=experiment_n_iter)
    except attractors.TimeoutError as e:
        model_bitchange_probability = numpy.inf
    print "time taken for model_bitchange_probability={:.2f} secs".format(time.time() - second_start)
    third_start = time.time()
    try:
        state_bitchange_probability = \
            attractors.find_state_bitchange_probability_for_different_attractors(graph, n_iter=experiment_n_iter)
    except attractors.TimeoutError as e:
        state_bitchange_probability = numpy.inf
    # print "time taken for state_bitchange_probability={:.2f} secs".format(time.time() - third_start)

    result = StabilityResult(graph_name=name, random_functions=False, random_edges=False,
                             n=len(graph.vertices), minimal_model_bitchange=minimal_model_bitchange,
                             model_bitchange_prob=model_bitchange_probability,
                             state_bitchange_prob=state_bitchange_probability)
    results.append(result)

    # "bitchange={},
    print "model_prob={:.2f}, state_prob={:.2f}, time_taken={:.2f} secs".format(
        minimal_model_bitchange,
        # model_bitchange_probability,
        state_bitchange_probability, time.time() - start
    )
for test in range(5000):
    print "test #{}".format(test + 1)
    for graph, name in zip(biological_graphs, biological_graph_names):
        start = time.time()
        graph_copy = graph.copy()

        randomize_functions, randomize_edges = random.choice([(False, True), (True, False), (True, True)])
        if randomize_functions:
            graph_copy.randomize_functions()
        if randomize_edges:
            graph_copy.randomize_edges()

        try:
            minimal_model_bitchange = \
                attractors.find_model_bitchange_for_new_attractor(graph_copy, max_len=experiment_max_len,
                                                                  verbose=False, use_dubrova=True,
                                                                  simplify_boolean=False)
        except attractors.TimeoutError as e:
            minimal_model_bitchange = numpy.inf
        # print "time taken for model_bitchange={:.2f} secs".format(time.time() - start)
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
                                                                                     n_iter=experiment_n_iter)
        except attractors.TimeoutError as e:
            state_bitchange_probability = numpy.inf
        # print "time taken for state_bitchange_probability={:.2f} secs".format(time.time() - third_start)

        result = StabilityResult(graph_name=name, random_functions=randomize_functions,
                                 random_edges=randomize_edges, n=len(graph_copy.vertices),
                                 minimal_model_bitchange=minimal_model_bitchange,
                                 model_bitchange_prob=model_bitchange_probability,
                                 state_bitchange_prob=state_bitchange_probability)
        # "bitchange={},
        results.append(result)
        print "model_prob={:.2f}, state_prob={:.2f}, time_taken={:.2f} secs".format(
            minimal_model_bitchange,
            # model_bitchange_probability,
            state_bitchange_probability, time.time() - start
        )

    # save on each iteration, why not
    with open('temp_stability_dict_with_model_bitchange_prob.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["graph_name", "random_functions", "random_edges", "n", "minimal_model_bitchange", "model_bitchange_prob",
                         "state_bitchange_prob"])
        for stability_result in results:
            writer.writerow(stability_result)

