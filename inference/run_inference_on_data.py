import os
from attractor_learning import graphs
import numpy as np
import shutil
import logging
import time
from inference import dummy_inference, binary_inference_ideas
from sklearn.model_selection import train_test_split
from validation import inference_scoring
import time
import functools

inference_method = dummy_inference.dummy_inference_method
# inference_method = binary_inference_ideas.infer_known_topology_general
# inference_method = binary_inference_ideas.infer_known_topology_symmetric
# inference_method = binary_inference_ideas.infer_unknown_topology_symmetric

data_dir = "../data/generated/different_data_generation_parameters_small/general_funcs_true_graph_no_sample_noise"
timestr = time.strftime("%Y%m%d-%H%M%S")
output_parent_dir = os.path.join("../inferred_models/different_data_generation_parameters_small", "{}_on_{}".format(
    inference_method.__name__, os.path.split(data_dir)[-1]))
train_size = 0.8
model_inference_timeout_secs = 5
allow_additional_edges = True
included_edges_relative_weight = -0.1
added_edges_relative_weight = -0.1


timestr = time.strftime("%Y%m%d-%H%M%S")
if os.path.exists(output_parent_dir):
    shutil.rmtree(output_parent_dir)
os.makedirs(output_parent_dir)

# concatenate log for model generation with the inference log
with open(os.path.join(output_parent_dir, "log.txt"), 'wb') as new_log_file:
    with open(os.path.join(data_dir, "log.txt"), 'rb') as old_log_file:
        shutil.copyfileobj(old_log_file, new_log_file)

logger = logging.getLogger()
logging.basicConfig(filename=os.path.join(output_parent_dir, "log.txt"),
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)


logger.info("inference_method_name={}".format(inference_method.__name__))
logger.info("inference_method_module={}".format(inference_method.__module__))
logger.info("train_size={}".format(train_size))
logger.info("model_inference_timeout_secs={}".format(model_inference_timeout_secs))
logger.info("allow_additional_edges={}".format(allow_additional_edges))
logger.info("included_edges_relative_weight={}".format(included_edges_relative_weight))
logger.info("added_edges_relative_weight={}".format(added_edges_relative_weight))


def main():

    paths = [(f.name, f.path) for f in os.scandir(data_dir) if f.is_dir()]
    if len(paths) == 0:
        raise ValueError("Did not find network directories in data_dir {}".format(data_dir))
    for network_name, network_path in paths:
        scaffold_network = graphs.Network.parse_cnet(os.path.join(network_path, "scaffold_network.cnet"))
        true_network = graphs.Network.parse_cnet(os.path.join(network_path, "true_network.cnet"))
        reference_matrices = np.load(os.path.join(network_path, "matrices.npz"))

        if train_size != 1:
            train_keys, test_keys = train_test_split(list(reference_matrices.keys()), train_size=train_size)
            reference_train = {i: reference_matrices[i] for i in train_keys}
            reference_test = {i: reference_matrices[i] for i in test_keys}
        else:
            reference_test = dict()
            reference_train = reference_matrices

        start = time.time()
        inferred_model = inference_method(reference_train.values(), scaffold_network,
                                          timeout_secs=model_inference_timeout_secs,
                                          log_file=os.path.join(output_parent_dir, "inference_log.txt"),
                                          allow_additional_edges=allow_additional_edges,
                                          included_edges_relative_weight=included_edges_relative_weight,
                                          added_edges_relative_weight=added_edges_relative_weight)
        time_taken = time.time() - start
        os.makedirs(os.path.join(output_parent_dir, network_name))
        inferred_model.export_to_cnet(os.path.join(output_parent_dir, network_name,
                                                   "inferred_network.cnet"))

        for group, data in zip(['train', 'test'], [reference_train, reference_test]):
            inferred_matrices = {i: inferred_model.next_states(data_matrix[0, :], data_matrix.shape[0] - 1)
                                 for i, data_matrix in data.items()}
            inferred_matrices = {i: np.array(mat) for (i, mat) in inferred_matrices.items()}
            assert len(inferred_matrices) == len(data)
            assert all(pred_mat.shape[0] == data[i].shape[0] for i, pred_mat in inferred_matrices.items())
            np.savez(os.path.join(output_parent_dir, network_name, "{}_matrices".format(group)), **inferred_matrices)

            ref_vec = np.concatenate([data[i][1:, ].flatten() for i in inferred_matrices.keys()])
            pred_vec = np.concatenate([inferred_matrices[i][1:, ].flatten() for i in inferred_matrices.keys()])
            timeseries_score = inference_scoring.sparse_accuracy_score(ref_vec, pred_vec)
            np.save(os.path.join(output_parent_dir, network_name, "timeseries_accuracy_score_{}".format(group)),
                    timeseries_score)

        np.save(os.path.join(output_parent_dir, network_name, "inference_time"), time_taken)

        y_true, y_pred = inference_scoring.models_to_edge_vectors(true_network, inferred_model, use_sparse=True)
        edge_score = inference_scoring.sparse_accuracy_score(y_true, y_pred)
        np.save(os.path.join(output_parent_dir, network_name, "edge_accuracy_score"), edge_score)


if __name__ == "__main__":
    main()
