import os
import shutil
from attractor_learning import graphs
from synthetic_data_generation.graph_generation.our_methods import generate_random_graphs, generate_scaffold_network
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


def main():
    p = configargparse.ArgParser(default_config_files=['./config.txt'])
    p.add_argument('-c', '--config', required=False, is_config_file=True, help='config file path to override defaults')
    p.add_argument('--use_random_network', required=False, type=bool)
    p.add_argument('--experiments_per_network', required=False, type=int, action='append')
    p.add_argument('--graphs_dir', required=False, type=str)
    p.add_argument('--data_output_parent_dir', required=False, type=str)
    p.add_argument('--timepoints_per_experiment', required=False, type=int, action='append')
    p.add_argument('--state_sample_type', required=False, type=str)
    p.add_argument('--frequency_handling', required=False, type=str)
    p.add_argument('--sample_to_model_freq_ratio', required=False, type=float)
    p.add_argument('--state_noise_chance', required=False, type=float, action='append')
    p.add_argument('--frequency_noise_std', required=False, type=float)
    p.add_argument('--random_networks_per_reference', required=False, type=int)
    p.add_argument('--mutate_input_nodes', required=False, type=bool)
    p.add_argument('--preserve_truth_ratio', required=False, type=bool)
    p.add_argument('--function_type_restriction', required=False, type=str, action='append')
    p.add_argument('--preserve_input_nodes_on_add', required=False, type=bool)
    p.add_argument('--scaffold_network_added_edge_fraction', required=False, type=float, action='append')
    p.add_argument('--scaffold_network_removed_edge_fraction', required=False, type=float, action='append')
    options = p.parse_args()

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
        options_combination['function_type_restriction'] = FunctionTypeRestriction[
            options_combination['function_type_restriction']]

        # need to represent the argument combination as a string to use in filename. Need to extract name from enums.
        comb_str = str({k: (v.name if isinstance(v, enum.Enum) else v) for k, v in options_combination.items()})
        comb_str = comb_str.translate(str.maketrans('', '', "'{}")).replace(": ", "=").replace(", ", "_")

        kwargs = options_combination | constant_options

        data_dir_path = os.path.join(kwargs['data_output_parent_dir'], comb_str)
        if os.path.exists(data_dir_path):
            shutil.rmtree(data_dir_path)
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
            reference_graphs = generate_random_graphs(kwargs['graphs_dir'], options=options)
        else:
            reference_graphs = []
            for graph_dir in os.listdir(kwargs['graphs_dir']):
                G = graphs.Network.parse_boolean_tables(os.path.join(kwargs['graphs_dir'], graph_dir))
                reference_graphs.append(G)

        for graph_index, reference_graph in enumerate(reference_graphs):
            graph_path = os.path.join(data_dir_path, "network_{}".format(graph_index))
            os.makedirs(graph_path)
            reference_graph.export_to_cnet(os.path.join(graph_path, "true_network.cnet"))
            random_scaffold = generate_scaffold_network(reference_graph, **kwargs)
            for vertex in random_scaffold.vertices:
                vertex.function = lambda *x: False  # so it can fit cnet format.
            random_scaffold.export_to_cnet(os.path.join(graph_path, "scaffold_network.cnet"))
            matrices = generate_experiments_data(reference_graph, kwargs['experiments_per_network'], **kwargs)
            named_matrices = dict((str(i), mat) for (i, mat) in enumerate(matrices))
            np.savez(os.path.join(graph_path, "matrices"), **named_matrices)


if __name__ == "__main__":
    main()
