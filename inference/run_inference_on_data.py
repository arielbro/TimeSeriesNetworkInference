import os
from attractor_learning import graphs
import numpy as np
import shutil
import logging
import time
import enum
from inference import dummy_inference, binary_inference_ideas
from inference.binary_inference_ideas import infer_known_topology_symmetric, infer_known_topology_general, infer_unknown_topology_symmetric
from sklearn.model_selection import train_test_split
from validation import inference_scoring
import time
import re
import itertools
import functools
import configargparse


def main():
    p = configargparse.ArgParser(default_config_files=['./config.txt'])
    p.add_argument('--config', is_config_file=True, help='config file path to override defaults')
    p.add_argument('--inference_method', required=False, type=str, action='append')
    p.add_argument('--data_parent_dir', required=False, type=str, action='append')
    p.add_argument('--train_size', required=False, type=float)
    p.add_argument('--model_inference_timeout_secs', required=False, type=float)
    p.add_argument('--allow_additional_edges', required=False, default=False, type=bool)
    p.add_argument('--included_edges_relative_weight', required=False, type=float)
    p.add_argument('--added_edges_relative_weight', required=False, type=float)
    options = p.parse_args()

    constant_options = {k: v for (k, v) in options._get_kwargs() if not isinstance(v, list)}
    variable_options = {k: v for (k, v) in options._get_kwargs() if isinstance(v, list)}
    options_combinations = (dict(zip(variable_options, x)) for x in itertools.product(*variable_options.values()))

    # run over different combinations of options as specified in the config
    for options_combination in options_combinations:

        # need to translate inference method to the enum
        if options_combination['inference_method'] == "random":
            inference_method = dummy_inference.dummy_inference_method
        elif options_combination['inference_method'] == "general":
            inference_method = infer_known_topology_general
        elif options_combination['inference_method'] == "symmetric":
            inference_method = infer_known_topology_symmetric
        elif options_combination['inference_method'] == "symmetric_topology":
            inference_method = infer_unknown_topology_symmetric
        else:
            raise ValueError("Unknown inference method {}".format(options_combination['inference_method']))
        options_combination['inference_method'] = inference_method

        # need to represent the argument combination as a string to use in filename. Need to extract name from enums
        # and from inference functions
        comb_str = str({k: (v.name if isinstance(v, enum.Enum) else v) for k, v in options_combinations.items()})
        comb_str = comb_str.translate(str.maketrans('', '', "'{}")).replace(": ", "=").replace(", ", "_")
        comb_str = re.sub('<function ([a-zA-Z_]+) at.*', '\\1', comb_str)

        # kwargs = options_combination | constant_options (works on python>=3.9)
        kwargs = options_combination.copy()
        kwargs.update(constant_options)

        for data_dir in os.scandir(kwargs['data_parent_dir']):
            if not os.path.isdir(data_dir):
                continue
            output_parent_dir = os.path.join("inferred_models", os.path.split(kwargs['data_parent_dir'])[-1],
                "{}_data_dir={}".format(comb_str, data_dir.name))
            if os.path.exists(output_parent_dir):
                shutil.rmtree(output_parent_dir, ignore_errors=True)
            os.makedirs(output_parent_dir)

            # concatenate log for model generation with the inference log
            with open(os.path.join(output_parent_dir, "log.txt"), 'wb') as new_log_file:
                with open(os.path.join(data_dir.path, "log.txt"), 'rb') as old_log_file:
                    shutil.copyfileobj(old_log_file, new_log_file)

            logger = logging.getLogger()
            logging.basicConfig(filename=os.path.join(output_parent_dir, "log.txt"),
                                filemode='a',
                                format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                                datefmt='%H:%M:%S',
                                level=logging.DEBUG)
            for var, val in kwargs.items():
                logger.info("{}={}".format(var, val))

            paths = [(f.name, f.path) for f in os.scandir(data_dir.path) if f.is_dir()]
            if len(paths) == 0:
                raise ValueError("Did not find network directories in data_dir {}".format(data_dir.name))
            for network_name, network_path in paths:
                scaffold_network = graphs.Network.parse_cnet(os.path.join(network_path, "scaffold_network.cnet"))
                true_network = graphs.Network.parse_cnet(os.path.join(network_path, "true_network.cnet"))
                reference_matrices = np.load(os.path.join(network_path, "matrices.npz"))

                if kwargs['train_size'] != 1:
                    train_keys, test_keys = train_test_split(list(reference_matrices.keys()),
                                                             train_size=kwargs['train_size'])
                    reference_train = {i: reference_matrices[i] for i in train_keys}
                    reference_test = {i: reference_matrices[i] for i in test_keys}
                else:
                    reference_test = dict()
                    reference_train = reference_matrices

                start = time.time()
                inferred_model = inference_method(reference_train.values(), scaffold_network,
                                                  timeout_secs=kwargs['model_inference_timeout_secs'],
                                                  log_file=os.path.join(output_parent_dir, "inference_log.txt"),
                                                  allow_additional_edges=kwargs['allow_additional_edges'],
                                                  included_edges_relative_weight=kwargs['included_edges_relative_weight'],
                                                  added_edges_relative_weight=kwargs['added_edges_relative_weight'])
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
