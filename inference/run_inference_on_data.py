import os
from attractor_learning import graphs
import numpy as np
import shutil
import logging
import time
import datetime
import enum
import traceback
from concurrent.futures import ProcessPoolExecutor
from inference import dummy_inference, binary_inference_ideas
from inference.binary_inference_ideas import infer_known_topology_symmetric, infer_known_topology_general, infer_unknown_topology_symmetric
from sklearn.model_selection import train_test_split
from validation import inference_scoring
import time
import re
import itertools
import functools
import configargparse


def process_network(network_name, network_path, output_parent_dir, kwargs):
    inference_method = kwargs['inference_method']
    network_out_dir = os.path.join(output_parent_dir, network_name)
    os.makedirs(network_out_dir, exist_ok=True)
    network_log_file = os.path.join(network_out_dir, "inference_log.txt")
    network_started_iso = datetime.datetime.now().isoformat(timespec='seconds')
    network_wall_t0 = time.time()
    status = 'ok'

    try:
        scaffold_network = graphs.Network.parse_cnet(os.path.join(network_path, "scaffold_network.cnet"))
        true_network = graphs.Network.parse_cnet(os.path.join(network_path, "true_network.cnet"))
        with np.load(os.path.join(network_path, "noisy_matrices.npz")) as reference_matrices, \
                np.load(os.path.join(network_path, "real_matrices.npz")) as hidden_matrices:
            if kwargs['train_size'] != 1:
                train_keys, test_keys = train_test_split(list(reference_matrices.keys()),
                                                         train_size=kwargs['train_size'])
                reference_train = {i: reference_matrices[i] for i in train_keys}
                real_train = {i: hidden_matrices[i] for i in train_keys}
                reference_test = {i: reference_matrices[i] for i in test_keys}
                real_test = {i: hidden_matrices[i] for i in test_keys}
            else:
                reference_test = dict()
                reference_train = {i: reference_matrices[i] for i in reference_matrices.keys()}
                real_test = dict()
                real_train = {i: hidden_matrices[i] for i in hidden_matrices.keys()}

        start = time.time()
        inferred_model = inference_method(reference_train.values(), scaffold_network,
                                          timeout_secs=kwargs['model_inference_timeout_secs'],
                                          log_file=network_log_file,
                                          allow_additional_edges=kwargs['allow_additional_edges'],
                                          included_edges_relative_weight=kwargs['included_edges_relative_weight'],
                                          added_edges_relative_weight=kwargs['added_edges_relative_weight'],
                                          allow_input_flips=kwargs['allow_input_flips'],
                                          flip_penalty=kwargs['flip_penalty'],
                                          no_anchoring=kwargs['no_anchoring'])
        time_taken = time.time() - start
        del scaffold_network
        inferred_model.export_to_cnet(os.path.join(network_out_dir, "inferred_network.cnet"))

        for group, reference_data, real_data in zip(['train', 'test'], [reference_train, reference_test], [real_train, real_test]):
            inferred_matrices = {i: np.asarray(inferred_model.next_states(data_matrix[0, :], data_matrix.shape[0] - 1))
                                 for i, data_matrix in reference_data.items()}
            assert len(inferred_matrices) == len(reference_data)
            assert all(pred_mat.shape[0] == reference_data[i].shape[0] for i, pred_mat in inferred_matrices.items())
            np.savez(os.path.join(network_out_dir, "{}_matrices".format(group)), **inferred_matrices)

            # compare against non-noisy data
            pred_vec = np.concatenate([inferred_matrices[i][1:, ].flatten() for i in inferred_matrices.keys()])
            ref_vec = np.concatenate([real_data[i][1:, ].flatten() for i in inferred_matrices.keys()])
            timeseries_score = inference_scoring.sparse_accuracy_score(ref_vec, pred_vec)
            np.save(os.path.join(network_out_dir, "timeseries_real_accuracy_score_{}".format(group)),
                    timeseries_score)
            del ref_vec
            # compare against given reference data
            ref_vec = np.concatenate([reference_data[i][1:, ].flatten() for i in inferred_matrices.keys()])
            timeseries_score = inference_scoring.sparse_accuracy_score(ref_vec, pred_vec)
            np.save(os.path.join(network_out_dir, "timeseries_reference_accuracy_score_{}".format(group)),
                    timeseries_score)
            del ref_vec, pred_vec, inferred_matrices

        del reference_train, real_train, reference_test, real_test

        np.save(os.path.join(network_out_dir, "inference_time"), time_taken)

        y_true, y_pred = inference_scoring.models_to_edge_vectors(true_network, inferred_model, use_sparse=True)
        edge_score = inference_scoring.sparse_jaccard_score(y_true, y_pred)
        np.save(os.path.join(network_out_dir, "edge_accuracy_score"), edge_score)
        del true_network, inferred_model, y_true, y_pred
    except Exception:
        status = 'error'
        with open(network_log_file, 'a', encoding='utf-8') as f:
            f.write("\nERROR processing network={} pid={}\n".format(network_name, os.getpid()))
            f.write(traceback.format_exc())

    return (network_name, network_log_file, network_started_iso,
            time.time() - network_wall_t0, status)


def append_network_log_to_master(master_log_path, network_name, network_log_path,
                                  started_iso, elapsed_s, status):
    with open(master_log_path, 'a', encoding='utf-8') as out:
        out.write("\n>>> begin network={name} started={started} elapsed={elapsed:.2f}s status={status}\n".format(
            name=network_name, started=started_iso, elapsed=elapsed_s, status=status))
        if os.path.exists(network_log_path) and os.path.getsize(network_log_path) > 0:
            with open(network_log_path, 'r', encoding='utf-8') as src:
                shutil.copyfileobj(src, out)
            out.write("\n")
        else:
            out.write("(no per-network log produced)\n")
        out.write("<<< end network={name}\n".format(name=network_name))


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
    p.add_argument('--allow_input_flips', required=False, default=False, type=bool)
    p.add_argument('--flip_penalty', required=False, default=1.0, type=float)
    p.add_argument('--no_anchoring', required=False, default=False, type=bool)
    p.add_argument('--n_processes', required=False, type=int, default=1)
    options = p.parse_args()

    constant_options = {k: v for (k, v) in options._get_kwargs() if not isinstance(v, list)}
    variable_options = {k: v for (k, v) in options._get_kwargs() if isinstance(v, list)}
    options_combinations = (dict(zip(variable_options, x)) for x in itertools.product(*variable_options.values()))

    # run over different combinations of options as specified in the config
    for options_combination in options_combinations:
        print("Current parameters: ", options_combination)

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
        comb_str = str({k: (v.name if isinstance(v, enum.Enum) else v) for k, v in options_combination.items()
                        if k != 'data_parent_dir'})
        comb_str = comb_str.translate(str.maketrans('', '', "'{}")).replace(": ", "=").replace(", ", "_")
        comb_str = re.sub('<function ([a-zA-Z_]+) at.*', '\\1', comb_str)

        # kwargs = options_combination | constant_options (works on python>=3.9)
        kwargs = options_combination.copy()
        kwargs.update(constant_options)

        for data_dir in os.scandir(kwargs['data_parent_dir']):
            print("Processing data in {}".format(data_dir.path))
            if not os.path.isdir(data_dir):
                continue
            data_dir_partial_path = os.path.join(*os.path.normpath(kwargs['data_parent_dir']).split(os.sep)[-2:])
            output_parent_dir = os.path.join("inferred_models", data_dir_partial_path,
                "{}-{}".format(comb_str, data_dir.name))
            if os.path.exists(output_parent_dir):
                shutil.rmtree(output_parent_dir, ignore_errors=True)
            os.makedirs(output_parent_dir)

            # concatenate log for model generation with the inference log
            # with open(os.path.join(output_parent_dir, "log.txt"), 'wb') as new_log_file:
            #     with open(os.path.join(data_dir.path, "log.txt"), 'rb') as old_log_file:
            #         shutil.copyfileobj(old_log_file, new_log_file)

            logger = logging.getLogger()
            for h in list(logger.handlers):
                logger.removeHandler(h)
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

            n_processes = max(1, int(kwargs.get('n_processes', 1) or 1))
            worker = functools.partial(process_network,
                                       output_parent_dir=output_parent_dir,
                                       kwargs=kwargs)
            if n_processes > 1:
                with ProcessPoolExecutor(max_workers=n_processes) as ex:
                    network_log_meta = list(ex.map(worker,
                                                   [name for name, _ in paths],
                                                   [path for _, path in paths]))
            else:
                network_log_meta = [worker(name, path) for name, path in paths]

            master_log_path = os.path.join(output_parent_dir, "log.txt")
            for h in logger.handlers:
                h.flush()
            for meta in network_log_meta:
                append_network_log_to_master(master_log_path, *meta)


if __name__ == "__main__":
    main()
