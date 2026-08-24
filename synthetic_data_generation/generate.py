import os
import shutil
from attractor_learning import graphs
from synthetic_data_generation.graph_generation.our_methods import generate_random_graphs, \
    generate_scaffold_network, parse_added_edge_fraction, randomize_reference_graph
from concurrent.futures import ProcessPoolExecutor
import random
from synthetic_data_generation.time_series_generation.our_methods import generate_experiments_data
from synthetic_data_generation.time_series_generation.our_methods import FrequencyHandling, StateSampleType
import time
import numpy as np
import logging
import enum
import configargparse
import json
import functools
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


# Files generate_data writes for one network. A directory holding all of them (non-empty) is finished and
# is skipped on a re-run.
NETWORK_OUTPUT_FILES = ("source_model.txt", "true_network.json", "scaffold_network.json",
                        "real_matrices.npz", "noisy_matrices.npz")

# Parameters recorded beside the generated data, so a re-run can tell whether it is continuing the same job.
# The config file's path and the output location are excluded: neither changes what gets generated, and
# requiring them to match would refuse a resume that is actually valid.
GENERATION_PARAMS_FILE = "generation_params.json"
# n_processes and seed are excluded as well: neither changes the distribution the data is drawn from,
# so a resume with a different process count or seed is continuing the same job.
_FINGERPRINT_EXCLUDED = {"config", "data_output_parent_dir", "regenerate", "n_processes", "seed"}


def network_data_is_complete(graph_path):
    """True iff graph_path already holds every (non-empty) file generation writes for one network, so it can
    be skipped. A network killed part-way through is missing at least one and is regenerated from scratch."""
    return all(os.path.exists(os.path.join(graph_path, name)) and
               os.path.getsize(os.path.join(graph_path, name)) > 0
               for name in NETWORK_OUTPUT_FILES)


def generation_fingerprint(kwargs):
    """JSON-serializable record of the parameters a data directory was generated with."""
    return {key: (value.name if isinstance(value, enum.Enum) else value)
            for key, value in sorted(kwargs.items()) if key not in _FINGERPRINT_EXCLUDED}


def prepare_data_dir(data_dir_path, kwargs, regenerate):
    """Make data_dir_path ready to be written into, and return whether finished networks in it may be kept.

    Resuming is only safe when the run is the same job: the directory name carries the swept parameters, but
    the constant ones (graphs_dir, state_sample_type, sample_to_model_freq_ratio, ...) are not in it, so a
    changed constant would otherwise mix data from two different configurations under one name. The
    parameters are therefore recorded in the directory and compared; anything that does not match, including
    a directory written before this file existed, is regenerated rather than continued.

    use_random_network never resumes, whatever the parameters say. There the reference graphs are randomized
    afresh on every run, so network_<index> is a different network each time: the networks kept from an
    earlier run would not be the ones this run produced, and the folder would end up holding graphs from two
    unrelated draws with nothing recording which came from where."""
    fingerprint = generation_fingerprint(kwargs)
    fingerprint_path = os.path.join(data_dir_path, GENERATION_PARAMS_FILE)
    if os.path.exists(data_dir_path) and not regenerate and not kwargs.get('use_random_network'):
        stored = None
        if os.path.exists(fingerprint_path):
            try:
                with open(fingerprint_path) as f:
                    stored = json.load(f)
            except ValueError:
                stored = None
        if stored == fingerprint:
            print("Continuing {}: finished networks in it will be kept".format(data_dir_path))
            return True
        reason = "it was generated with different parameters" if stored is not None else \
            "it holds no record of the parameters it was generated with"
        print("Regenerating {} from scratch: {}".format(data_dir_path, reason))
    if os.path.exists(data_dir_path):
        if regenerate:
            print("Regenerating {} from scratch: --regenerate was given".format(data_dir_path))
        elif kwargs.get('use_random_network'):
            print("Regenerating {} from scratch: use_random_network randomizes the reference graphs per "
                  "run, so networks from an earlier run cannot be continued".format(data_dir_path))
        shutil.rmtree(data_dir_path, ignore_errors=True)
    os.makedirs(data_dir_path)
    with open(fingerprint_path, 'w') as f:
        json.dump(fingerprint, f, indent=1)
    return False


def generate_one_network(task, kwargs):
    """Generate one network's data. Runs in a worker process, so it takes a description of the reference
    model rather than the model itself: a Network whose functions are truth tables holds a lambdified
    function and cannot be pickled. Parsing here also spreads the parse cost over the pool.

    The seed is per network and set explicitly. Worker processes started by fork inherit the parent's RNG
    state, so without this every worker would draw the same scaffold and the same starting states."""
    graph_path, model_dir, source_label, randomize, seed = task
    random.seed(seed)
    np.random.seed(seed % (2 ** 32))

    reference_graph = graphs.Network.parse_model_dir(model_dir)
    if randomize:
        reference_graph = randomize_reference_graph(reference_graph, **kwargs)

    os.makedirs(graph_path, exist_ok=True)
    # the index is positional, so record which model it came from - otherwise a folder of named
    # models (the toy families, say) is unreadable downstream
    with open(os.path.join(graph_path, "source_model.txt"), 'w') as source_file:
        source_file.write("{}\n".format(source_label))
    reference_graph.save(os.path.join(graph_path, "true_network.json"))
    random_scaffold = generate_scaffold_network(reference_graph, **kwargs)
    random_scaffold.save(os.path.join(graph_path, "scaffold_network.json"))

    named_real_matrices, named_noisy_matrices = dict(), dict()
    for i, (real_mat, noisy_mat) in enumerate(generate_experiments_data(reference_graph, **kwargs)):
        named_real_matrices[str(i)] = real_mat
        named_noisy_matrices[str(i)] = noisy_mat
    np.savez(os.path.join(graph_path, "real_matrices"), **named_real_matrices)
    np.savez(os.path.join(graph_path, "noisy_matrices"), **named_noisy_matrices)
    return graph_path


def reference_models(kwargs):
    """(relative path, model directory, randomize?) for every reference this run generates data from, in a
    fixed order, without parsing any of them - the index a network gets is positional, so the order has to
    be deterministic, but the parsing belongs in the workers.

    With use_random_network each source model contributes random_networks_per_reference randomizations, all
    in one group (a randomized model has no source directory to be named after)."""
    if kwargs['use_random_network']:
        found = []
        for name in sorted(os.listdir(kwargs['graphs_dir'])):
            model_dir = os.path.join(kwargs['graphs_dir'], name)
            if not os.path.isdir(model_dir):
                continue     # a graphs_dir may hold loose files next to the model directories
            if (kwargs['max_graph_size'] is not None) and \
                    (graphs.Network.model_dir_size(model_dir) > kwargs['max_graph_size']):
                continue
            found.extend((None, model_dir, True) for _ in range(kwargs['random_networks_per_reference']))
        return found

    found = []
    for relative_path, model_dir in find_model_dirs(kwargs['graphs_dir']):
        # max_graph_size bounds the number of nodes; None (the default, i.e. the option left out of
        # the config) means no filtering. Probed from the model's index/header so an oversized model
        # is never parsed - parse_boolean_tables costs 2**in-degree per node to build its functions.
        if (kwargs['max_graph_size'] is not None) and \
                (graphs.Network.model_dir_size(model_dir) > kwargs['max_graph_size']):
            continue
        found.append((relative_path, model_dir, False))
    return found


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
    # by default a re-run continues a data directory generated with the same parameters, skipping the
    # networks already finished in it; this forces the directory to be rebuilt from scratch instead
    p.add_argument('--regenerate', required=False, default=False, type=parse_bool_option)
    p.add_argument('--n_processes', required=False, type=int, default=1,
                   help='networks generated in parallel; each is independent, and the cost is dominated by the per-network time-series generation')
    p.add_argument('--seed', required=False, type=int, default=None,
                   help='base seed for the per-network seeds; omit for a different draw each run')
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
        resuming = prepare_data_dir(data_dir_path, kwargs, kwargs['regenerate'])

        logger = logging.getLogger()
        logging.basicConfig(filename=os.path.join(kwargs['data_output_parent_dir'], comb_str, "log.txt"),
                            filemode='a',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)
        for var, val in kwargs.items():
            logger.info("{}={}".format(var, val))

        # network_<index> numbering restarts inside each group, so a grouped graphs_dir keeps its shape
        # here: data_dir_path/<group>/network_<index>. Flat graphs_dirs (group "") are unchanged.
        references = reference_models(kwargs)
        seed_source = random.Random(kwargs['seed'])
        index_in_group = dict()
        n_skipped = 0
        tasks = []
        for relative_path, model_dir, randomize in references:
            group = os.path.dirname(relative_path) if relative_path else ""
            graph_index = index_in_group.get(group, 0)
            index_in_group[group] = graph_index + 1
            graph_path = os.path.join(data_dir_path, group, "network_{}".format(graph_index))
            seed = seed_source.randrange(2 ** 31)   # drawn for every network, so skipping does not shift it
            # the index is positional, so it has to be advanced for every reference graph whether or not
            # this one is regenerated - hence the skip here rather than earlier
            if resuming and network_data_is_complete(graph_path):
                logger.info("skipping {}: already generated".format(graph_path))
                n_skipped += 1
                continue
            # '/' regardless of platform: source_model.txt is written on whichever machine generated the
            # data and read on whichever machine analyses it
            source_label = relative_path.replace(os.sep, "/") if relative_path \
                else "randomized_network_{}".format(graph_index)
            tasks.append((graph_path, model_dir, source_label, randomize, seed))

        n_processes = max(1, int(kwargs.get('n_processes', 1) or 1))
        if tasks and n_processes > 1:
            print("Generating {} networks over {} processes".format(len(tasks), n_processes))
            with ProcessPoolExecutor(max_workers=n_processes) as executor:
                list(executor.map(functools.partial(generate_one_network, kwargs=kwargs), tasks))
        else:
            for task in tasks:
                generate_one_network(task, kwargs)

        if n_skipped:
            print("Kept {} of {} networks already generated in {}".format(
                n_skipped, len(references), data_dir_path))


if __name__ == "__main__":
    main()
