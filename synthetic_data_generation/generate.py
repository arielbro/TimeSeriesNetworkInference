import os
import shutil
from attractor_learning import graphs
from graph_generation.our_methods import generate_random_graphs
from time_series_generation.our_methods import generate_one_experiment_data
import time
import numpy as np
import logging

experiments_per_network = 40
data_output_parent_dir = "../data/generated"

timestr = time.strftime("%Y%m%d-%H%M%S")
data_dir_path = os.path.join(data_output_parent_dir, timestr)
if os.path.exists(data_dir_path):
    shutil.rmtree(data_dir_path)
os.makedirs(data_dir_path)

logger = logging.getLogger(__name__)
logging.basicConfig(filename=os.path.join(data_dir_path, "log"),
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)


logger.info("experiments_per_network={}".format(experiments_per_network))


def main():

    data = []

    random_graphs = generate_random_graphs()

    for graph_index, random_graph in enumerate(random_graphs):
        graph_path = os.path.join(data_dir_path, "network_{}".format(graph_index))
        os.makedirs(graph_path)
        for experiment_index in range(experiments_per_network):
            matrix = generate_one_experiment_data(random_graph)
            matrix = np.array(matrix, dtype=bool)
            np.savetxt(os.path.join(graph_path, str(experiment_index)), matrix, delimiter=",",
                       fmt="%i")

    for i, matrix in enumerate(data):
        matrix = np.array(matrix)
        np.savetxt(os.path.join(data_dir_path, i), matrix, delimeter=",")


if __name__ == "__main__":
    main()
