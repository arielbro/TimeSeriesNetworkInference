import time
import graphs
import attractors
from collections import namedtuple
import csv

StabilityResult = namedtuple("StabilityResult", "graph_name is_random n minimal_model_bitchange "
                                                "model_bitchange_prob state_bitchange_prob")

results = []
biological_graphs = [
    graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
                              "BNS_Dubrova_2011\\MAPK_large2.cnet"),
    graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
                              "BNS_Dubrova_2011\\arabidopsis.cnet.txt"),
    graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
                              "BNS_Dubrova_2011\\EGFR_man.cnet"),
    graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
                              "BNS_Dubrova_2011\\MAPK_large.cnet"),
    graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
                              "BNS_Dubrova_2011\\thelper.cnet.txt"),
    graphs.Network.parse_cnet("C:\\Users\\ariel\\Downloads\\Attractors - for Ariel\\Attractors - for Ariel\\"
                              "BNS_Dubrova_2011\\tcr.cnet"),
]
biological_graph_names = ["MAPK_large2", "arabidopsis", "EGFR_man", "MAPK_large", "thelper", "tcr"]

for graph, name in zip(biological_graphs, biological_graph_names):
    minimal_model_bitchange = \
        attractors.find_model_bitchange_for_new_attractor(graph, max_len=15, verbose=False, use_dubrova=True)
    model_bitchange_probability = \
        attractors.find_model_bitchange_probability_for_different_attractors(graph, n_iter=10)
    state_bitchange_probability = \
        attractors.find_state_bitchange_probability_for_different_attractors(graph, n_iter=10)
    result = StabilityResult(graph_name=name, is_random=False, n=len(graph.vertices),
                             minimal_model_bitchange=minimal_model_bitchange,
                             model_bitchange_prob=model_bitchange_probability,
                             state_bitchange_prob=state_bitchange_probability)
    results.append(result)
    print "bitchange={}, model_prob={}, state_prob={}, time_taken={:.2f} secs".format(
        minimal_model_bitchange, model_bitchange_probability, state_bitchange_probability, time.time() - start
    )

for test in range(5000):
    for graph, name in zip(biological_graphs, biological_graph_names):
        start = time.time()
        graph.randomize_functions()
        minimal_model_bitchange = \
            attractors.find_model_bitchange_for_new_attractor(graph, max_len=15, verbose=False, use_dubrova=True)
        model_bitchange_probability = \
            attractors.find_model_bitchange_probability_for_different_attractors(graph, n_iter=10)
        state_bitchange_probability = \
            attractors.find_state_bitchange_probability_for_different_attractors(graph, n_iter=10)
        result = StabilityResult(graph_name=name, is_random=True, n=len(graph.vertices),
                                 minimal_model_bitchange=minimal_model_bitchange,
                                 model_bitchange_prob=model_bitchange_probability,
                                 state_bitchange_prob=state_bitchange_probability)
        print "bitchange={}, model_prob={}, state_prob={}, time_taken={:.2f} secs".format(
            minimal_model_bitchange, model_bitchange_probability, state_bitchange_probability, time.time() - start
        )

    # save on each iteration, why not
    with open('temp_stability_dict.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["graph_name", "is_random", "n", "minimal_model_bitchange", "model_bitchange_prob",
                         "state_bitchange_prob"])
        for stability_result in results:
            writer.writerow(stability_result)

