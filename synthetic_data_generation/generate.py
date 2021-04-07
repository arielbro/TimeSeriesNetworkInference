import os
import shutil
from attractor_learning import graphs
from synthetic_data_generation.graph_generation.our_methods import generate_random_graphs, generate_scaffold_network
from synthetic_data_generation.graph_generation.our_methods import log_params as graph_generation_log_params
from synthetic_data_generation.time_series_generation.our_methods import generate_experiments_data
from synthetic_data_generation.time_series_generation.our_methods import log_params as time_series_generation_log_params
import time
import numpy as np
import logging

experiments_per_network = 30
data_output_parent_dir = "../data/generated"

name_prefix = "general_funcs_false_graph_no_sample_noise"
# timestr = time.strftime("%Y%m%d-%H%M%S")
data_dir_path = os.path.join(data_output_parent_dir, name_prefix)
if os.path.exists(data_dir_path):
    shutil.rmtree(data_dir_path)
os.makedirs(data_dir_path)

logger = logging.getLogger()
logging.basicConfig(filename=os.path.join(data_dir_path, "log.txt"),
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)


logger.info("experiments_per_network={}".format(experiments_per_network))
graph_generation_log_params()
time_series_generation_log_params()


def main():

    random_graphs = generate_random_graphs(logger=logger)

    for graph_index, random_graph in enumerate(random_graphs):
        graph_path = os.path.join(data_dir_path, "network_{}".format(graph_index))
        os.makedirs(graph_path)
        random_graph.export_to_cnet(os.path.join(graph_path, "true_network.cnet"))
        random_scaffold = generate_scaffold_network(random_graph, logger=logger)
        for vertex in random_scaffold.vertices:
            vertex.function = lambda *x: False  # so it can fit cnet format.
        random_scaffold.export_to_cnet(os.path.join(graph_path, "scaffold_network.cnet"))
        matrices = generate_experiments_data(random_graph, experiments_per_network, logger=logger)
        named_matrices = dict((str(i), mat) for (i, mat) in enumerate(matrices))
        np.savez(os.path.join(graph_path, "matrices"), **named_matrices)


if __name__ == "__main__":
    main()
