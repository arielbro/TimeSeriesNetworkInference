import os
import shutil
from attractor_learning import graphs
from synthetic_data_generation.graph_generation.our_methods import generate_random_graphs, \
    generate_scaffold_network, parse_added_edge_fraction
from synthetic_data_generation.time_series_generation.our_methods import generate_experiments_data
from synthetic_data_generation.time_series_generation.our_methods import FrequencyHandling, StateSampleType
import time
import numpy as np
import logging
import enum
import configargparse
from attractor_learning.graphs import FunctionTypeRestriction
import itertools
import re


def parse_bool_option(value):
    """configargparse `type` for a boolean given as a value rather than a bare flag (`only_attractors = False`
    in a config file). `type=bool` cannot be used for this: it just calls bool() on the string, so every
    non-empty value - "False" and "0" included - comes out True. Twin of the helper of the same name in
    inference/run_inference_on_data.py (kept local so the data generator doesn't import the inference stack)."""
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in ('true', 'yes', 'y', '1'):
        return True
    if normalized in ('false', 'no', 'n', '0'):
        return False
    raise ValueError("expected a boolean (true/false), got {!r}".format(value))


def find_model_dirs(graphs_dir):
    """(relative path, absolute path) for every model directory under graphs_dir, at whatever depth it
    sits: a graphs_dir may group its models in subfolders (and hold loose files like a README next to
    them), and that grouping is mirrored into the generated data tree. A directory that holds a model is
    not descended into.

    Sorted, because the generated folders are named network_<index> by position - os.listdir order is
    arbitrary and differs between platforms, so an unsorted walk would map the same index to different
    models on different machines."""
    found = []
    for root, dirs, _ in os.walk(graphs_dir):
        dirs.sort()
        if os.path.normpath(root) != os.path.normpath(graphs_dir) and graphs.Network.is_model_dir(root):
            found.append((os.path.relpath(root, graphs_dir), root))
            dirs[:] = []
    return sorted(found)


def main():
    p = configargparse.ArgParser(default_config_files=['./config.txt'])
    p.add_argument('-c', '--config', required=False, is_config_file=True, help='config file path to override defaults')
    p.add_argument('--use_random_network', required=False, default=False,  type=parse_bool_option)
    p.add_argument('--experiments_per_network', required=False, type=int, action='append')
    p.add_argument('--graphs_dir', required=False, type=str)
    p.add_argument('--max_graph_size', required=False, default=None, type=int, action='append',
                   help='skip source models with more than this many nodes; leave the option out of the '
                        'config (the default) to load every model regardless of size')
    p.add_argument('--data_output_parent_dir', required=False, type=str)
    p.add_argument('--timepoints_per_experiment', required=False, type=int, action='append')
    p.add_argument('--state_sample_type', required=True, type=str)
    p.add_argument('--frequency_handling', required=False, type=str)
    p.add_argument('--sample_to_model_freq_ratio', required=False, type=float)
    p.add_argument('--state_noise_chance', required=False, type=float, action='append')
    p.add_argument('--frequency_noise_std', required=False, type=float)
    p.add_argument('--random_networks_per_reference', required=False, type=int)
    p.add_argument('--mutate_input_nodes', required=False, default=False,  type=parse_bool_option)
    p.add_argument('--preserve_truth_ratio', required=False, default=False,  type=parse_bool_option)
    p.add_argument('--function_type_restriction', required=False, type=str, action='append')
    p.add_argument('--preserve_input_nodes_on_add', required=False, default=False,  type=parse_bool_option)
    p.add_argument('--scaffold_network_added_edge_fraction', required=False,
                   type=parse_added_edge_fraction, action='append',
                   help='fraction of the reference graph\'s edges to add to the scaffold, or "all" for a '
                        'scaffold holding every possible edge (which ignores '
                        '--scaffold_network_removed_edge_fraction)')
    p.add_argument('--scaffold_network_removed_edge_fraction', required=False, type=float, action='append')
    # appendable so a config can grid-search it (`only_attractors = [False, True]`), which also puts it in
    # comb_str and hence in the data directory name. The default is applied after parsing, since argparse
    # appends onto a list default instead of replacing it.
    p.add_argument('--only_attractors', required=False, default=None, type=parse_bool_option,
                   action='append')
    p.add_argument('--attractor_estimation_n_walks', required=False, type=int)
    options = p.parse_args()
    if options.only_attractors is None:
        options.only_attractors = [False]  # kept a list so it stays a (single-valued) varying option

    print(options)

    # Abbreviations for parameter names in the data directory name; without them a folder covering several
    # varying parameters gets unwieldy (and pushes toward path-length limits).
    name_replacements = {
        "experiments_per_network": "exppernet",
        "scaffold_network_added_edge_fraction": "scaffaddfrac",
        "timepoints_per_experiment": "tpperexp",
        "function_type_restriction": "fncrestriction",
        "scaffold_network_removed_edge_fraction": "scaffremovedfrac",
        "state_noise_chance": "statenoise",
        "max_graph_size": "maxsize",
        "only_attractors": "onlyattr",
    }

    constant_options = {k: v for (k, v) in options._get_kwargs() if not isinstance(v, list)}
    variable_options = {k: v for (k, v) in options._get_kwargs()  if isinstance(v, list)}
    options_combinations = (dict(zip(variable_options, x)) for x in itertools.product(*variable_options.values()))

    constant_options['state_sample_type'] = StateSampleType[
        constant_options['state_sample_type']]
    constant_options['frequency_handling'] = FrequencyHandling[
        constant_options['frequency_handling']]

    # run over different combinations of options as specified in the config
    for options_combination in options_combinations:

        print(options_combinations)
        if 'function_type_restriction' in options_combination:
            options_combination['function_type_restriction'] = FunctionTypeRestriction[
                options_combination['function_type_restriction']]

        # need to represent the argument combination as a string to use in filename. Need to extract name from enums.
        comb_str = str({k: (v.name if isinstance(v, enum.Enum) else v) for k, v in options_combination.items()})
        comb_str = comb_str.translate(str.maketrans('', '', "'{}")).replace(": ", "=").replace(", ", "-")

        # also replace the variable names with shorter versions using name_replacements
        for long_name, short_name in name_replacements.items():
            comb_str = comb_str.replace(long_name, short_name)


        # kwargs = options_combination | constant_options (works on python>=3.9)
        kwargs = options_combination.copy()
        kwargs.update(constant_options)

        data_dir_path = os.path.join(kwargs['data_output_parent_dir'], comb_str)
        if os.path.exists(data_dir_path):
            shutil.rmtree(data_dir_path, ignore_errors=True)
        os.makedirs(data_dir_path)

        logger = logging.getLogger()
        logging.basicConfig(filename=os.path.join(kwargs['data_output_parent_dir'], comb_str, "log.txt"),
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)
        for var, val in kwargs.items():
            logger.info("{}={}".format(var, val))

        if kwargs['use_random_network']:
            # randomized models have no source directory to name them after, so they all share one group
            reference_graphs = [(None, G) for G in generate_random_graphs(**kwargs)]
        else:
            reference_graphs = []
            for relative_path, graph_path in find_model_dirs(kwargs['graphs_dir']):
                # max_graph_size bounds the number of nodes; None (the default, i.e. the option left out of
                # the config) means no filtering. Probed from the model's index/header so an oversized model
                # is never parsed - parse_boolean_tables costs 2**in-degree per node to build its functions.
                if (kwargs['max_graph_size'] is not None) and \
                        (graphs.Network.model_dir_size(graph_path) > kwargs['max_graph_size']):
                    continue
                reference_graphs.append((relative_path, graphs.Network.parse_model_dir(graph_path)))

        # network_<index> numbering restarts inside each group, so a grouped graphs_dir keeps its shape
        # here: data_dir_path/<group>/network_<index>. Flat graphs_dirs (group "") are unchanged.
        index_in_group = dict()
        for relative_path, reference_graph in reference_graphs:
            group = os.path.dirname(relative_path) if relative_path else ""
            graph_index = index_in_group.get(group, 0)
            index_in_group[group] = graph_index + 1
            graph_path = os.path.join(data_dir_path, group, "network_{}".format(graph_index))
            os.makedirs(graph_path)
            # the index is positional, so record which model it came from - otherwise a folder of named
            # models (the toy families, say) is unreadable downstream
            with open(os.path.join(graph_path, "source_model.txt"), 'w') as source_file:
                # '/' regardless of platform: this file is written on whichever machine generated the data
                # and read on whichever machine analyses it
                source_file.write("{}\n".format(relative_path.replace(os.sep, "/") if relative_path
                                                else "randomized_network_{}".format(graph_index)))
            reference_graph.save(os.path.join(graph_path, "true_network.json"))
            random_scaffold = generate_scaffold_network(reference_graph, **kwargs)
            random_scaffold.save(os.path.join(graph_path, "scaffold_network.json"))
            matrices = generate_experiments_data(reference_graph, **kwargs)
            named_real_matrices = dict()
            named_noisy_matrices = dict()
            for (i, (real_mat, noisy_mat)) in enumerate(matrices):
                named_real_matrices[str(i)] = real_mat
                named_noisy_matrices[str(i)] = noisy_mat

            np.savez(os.path.join(graph_path, "real_matrices"), **named_real_matrices)
            np.savez(os.path.join(graph_path, "noisy_matrices"), **named_noisy_matrices)


if __name__ == "__main__":
    main()
