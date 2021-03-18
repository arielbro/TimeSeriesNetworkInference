import os
from attractor_learning import graphs
import numpy as np
import shutil
import logging
import time
from inference import dummy_inference

data_dir = "../data/generated/20210317-222947"
output_parent_dir = os.path.join("../inferred_models", os.path.split(data_dir)[-1])
inference_method = dummy_inference.dummy_inference_method

timestr = time.strftime("%Y%m%d-%H%M%S")
if os.path.exists(output_parent_dir):
    shutil.rmtree(output_parent_dir)
os.makedirs(output_parent_dir)

# concatenate log for model generation with the inference log
with open(os.path.join(output_parent_dir, "log"), 'wb') as new_log_file:
    with open(os.path.join(data_dir, "log"), 'rb') as old_log_file:
        shutil.copyfileobj(old_log_file, new_log_file)

logger = logging.getLogger()
logging.basicConfig(filename=os.path.join(output_parent_dir, "log"),
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)


logger.info("inference_method_name={}".format(inference_method.__name__))
logger.info("inference_method_module={}".format(inference_method.__module__))


def main():

    paths = [(f.name, f.path) for f in os.scandir(data_dir) if f.is_dir()]
    for network_name, network_path in paths:
        scaffold_network = graphs.Network.parse_cnet(os.path.join(network_path, "scaffold_network.cnet"))
        reference_matrices = np.load(os.path.join(network_path, "matrices.npz"))

        inferred_model = inference_method(reference_matrices.values(), scaffold_network)
        inferred_matrices = {i: inferred_model.next_states(data_matrix[0, :], data_matrix.shape[0] - 1)
                             for i, data_matrix in reference_matrices.items()}
        inferred_matrices = {i: np.array(mat) for (i, mat) in inferred_matrices.items()}
        os.makedirs(os.path.join(output_parent_dir, network_name))
        inferred_model.export_to_cnet(os.path.join(output_parent_dir, network_name,
                                                   "inferred_network.cnet"))

        np.savez(os.path.join(output_parent_dir, network_name, "matrices"), **inferred_matrices)


if __name__ == "__main__":
    main()
