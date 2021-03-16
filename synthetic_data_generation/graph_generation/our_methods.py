import os
import shutil
from attractor_learning import graphs
from synthetic_data_generation.time_series_generation.our_methods import \
    generate_one_experiment_data, StateSampleType, FrequencyHandling

import logging
logger = logging.getLogger(__name__)

random_networks_per_reference = 5
graphs_dir = "../attractor_learning/cellcollective_models"
mutate_input_nodes = True
preserve_truth_ratio = True


logger.info("random_networks_per_reference={}".format(random_networks_per_reference))
logger.info("graphs_dir={}".format(graphs_dir))
logger.info("mutate_input_nodes={}".format(mutate_input_nodes))
logger.info("preserve_truth_ratio={}".format(preserve_truth_ratio))

reference_graphs = []
for graph_dir in os.listdir(graphs_dir):
    G = graphs.Network.parse_boolean_tables(os.path.join(graphs_dir, graph_dir))
    reference_graphs.append(G)


def generate_random_graphs():

    for G in reference_graphs:
        for _ in range(random_networks_per_reference):
            random_graph = G.copy()
            random_graph.randomize_edges_by_switching(n_attempts=1000)
            random_graph.randomize_functions(mutate_input_nodes=mutate_input_nodes,
                                                            preserve_truth_ratio=preserve_truth_ratio)
            assert random_graph != G
            yield random_graph

