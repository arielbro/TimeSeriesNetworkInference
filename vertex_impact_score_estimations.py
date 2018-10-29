import cProfile
import random
import time
import graphs
import attractors
from shutil import copyfile
# from pebble import ProcessPool
from concurrent.futures import TimeoutError
from collections import namedtuple
import csv
import stochastic
import multiprocessing
import itertools
import numpy
import os
import sys
import signal, os, errno
from functools import wraps
import stat
import platform

stochastic_n_iter = 60
parallel = True
n_processes = 40
timeout_seconds = int(4 * 60 * 60)
graph_parent_dir = "cellcollective_models"
optimization_max_len = 1
optimization_max_num = 30
optimization_max_transient_len = 20
only_random = False
attractor_estimation_n_iter = 100
filter_out_timed_out_graphs = False
graph_size_filter = 35
n_attractors_filter = 50
queue_all_tasks = True
n_tests = 10

VertexImpactResult = namedtuple("VertexImpactResult", "graph_name is_random size "
                                                      "maximal_change_bits n_inputs normalized_n_inputs "
                                                      "max_degree mean_degree "
                                                      "num_attractors mean_attractor_length "
                                                      "median_attractor_length std_attractor_length "
                                                      "max_attractor_length std_attractor_basin "
                                                      "median_attractor_basin max_attractor_basin "
                                                      "optimization_model_impact_scores "
                                                      "stochastic_model_impact_scores "
                                                      "optimization_model_addition_impact_scores "
                                                      "stochastic_model_addition_impact_scores "
                                                      "optimization_state_impact_scores "
                                                      "stochastic_state_impact_scores "
                                                      "optimization_model_time "
                                                      "stochastic_model_time "
                                                      "optimization_model_addition_time "
                                                      "stochastic_model_addition_time "
                                                      "optimization_state_time "
                                                      "stochastic_state_time"
                                )


def timeout(seconds, error_message=os.strerror(errno.ETIMEDOUT)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.setitimer(signal.ITIMER_REAL, seconds)
            try:
                ret = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return ret

        return wraps(func)(wrapper)

    return decorator


def one_graph_impact_score_estimation_wrapper(args):
    try:
        return one_graph_impact_score_estimation(*args)
    except TimeoutError as e:
        print "warning - timeout on vertex impact score estimation after {} seconds.".format(timeout_seconds)
        return None
    except attractors.TimeoutError as e:
        print "Breaking impact score estimation (guorbi timeout)". \
            format(timeout_seconds)
        return None

@timeout(timeout_seconds)
def one_graph_impact_score_estimation(graph, name, is_biological, graph_name_to_attributes):
    # copy Dubrova's executable, to prevent mysterious errors when running multiple processes calling it.
    # process_specific_dubrova_path = "temp_{}_{}".format(attractors.dubrova_path, os.getpid())
    # copyfile(attractors.dubrova_path, process_specific_dubrova_path)
    # os.chmod(process_specific_dubrova_path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
    # attractors.dubrova_path = process_specific_dubrova_path
    # TODO: tidy up after done

    start = time.time()
    graph_copy = graph.copy()

    if is_biological:
        randomize_functions, randomize_edges = False, False
    else:
        # randomize_functions, randomize_edges = random.choice([(False, True), (True, False), (True, True)])
        randomize_functions, randomize_edges = True, True

    if randomize_functions:
        graph_copy.randomize_functions(preserve_truth_ratio=True)
    if randomize_edges:
        # graph_copy.randomize_incoming_edges()
        start = time.time()
        graph_copy.randomize_edges_by_switching()
    if randomize_edges or randomize_functions:
        print "time taken for graph randomization={:.2f} secs".format(time.time() - start)
    start = time.time()

    attractor_basin_tuples = stochastic.estimate_attractors(graph_copy, n_walks=1000, max_walk_len=100,
                                                            with_basins=True)
    current_attractors = [pair[0] for pair in attractor_basin_tuples]
    basin_sizes = [pair[1] for pair in attractor_basin_tuples]
    sum_sizes = sum(basin_sizes)
    basin_sizes = [size / float(sum_sizes) for size in basin_sizes]

    num_attractors = len(current_attractors)
    mean_attractor_length = numpy.mean([len(a) for a in current_attractors])
    median_attractor_length = numpy.median([len(a) for a in current_attractors])
    std_attractor_length = numpy.std([len(a) for a in current_attractors])
    max_attractor_length = numpy.max([len(a) for a in current_attractors])
    std_attractor_basin = numpy.std(basin_sizes)
    median_attractor_basin = numpy.median(basin_sizes)
    max_attractor_basin = numpy.max(basin_sizes)
    print "time taken for attractor estimation={:.2f} secs".format(time.time() - start)

    res_start = time.time()
    stochastic_model_impact_scores = attractors. \
        stochastic_vertex_model_impact_scores(graph_copy, current_attractors,
                                              use_dubrova=False,
                                              n_iter=stochastic_n_iter,
                                              impact_type=attractors.ImpactType.Invalidation,
                                              bits_of_change=1,
                                              attractor_estimation_n_iter=attractor_estimation_n_iter,
                                              relative_attractor_basin_sizes=basin_sizes)
    stochastic_model_time = time.time() - res_start

    res_start = time.time()
    optimization_model_impact_scores = attractors. \
        vertex_model_impact_scores(graph_copy, current_attractors=current_attractors,
                                   max_len=50, max_num=50, verbose=False,
                                   impact_types=attractors.ImpactType.Invalidation,
                                   normalize_addition_scores=True,
                                   relative_attractor_basin_sizes=basin_sizes,
                                   maximal_bits_of_change=1)
    optimization_model_time = time.time() - res_start

    res_start = time.time()
    stochastic_model_addition_impact_scores = attractors. \
        stochastic_vertex_model_impact_scores(graph_copy, current_attractors,
                                              use_dubrova=False,
                                              attractor_estimation_n_iter=attractor_estimation_n_iter,
                                              n_iter=stochastic_n_iter,
                                              bits_of_change=1,
                                              impact_type=attractors.ImpactType.Addition,
                                              relative_attractor_basin_sizes=basin_sizes)
    stochastic_model_addition_time = time.time() - res_start

    res_start = time.time()
    optimization_model_addition_impact_scores = attractors. \
        vertex_model_impact_scores(graph_copy, current_attractors=current_attractors,
                                   max_len=optimization_max_len, max_num=optimization_max_num, verbose=False,
                                   impact_types=attractors.ImpactType.Addition,
                                   normalize_addition_scores=True,
                                   relative_attractor_basin_sizes=basin_sizes,
                                   maximal_bits_of_change=1)
    optimization_model_addition_time = time.time() - res_start

    res_start = time.time()
    stochastic_state_impact_scores = attractors. \
        stochastic_vertex_state_impact_scores(graph_copy, n_iter=stochastic_n_iter)
    stochastic_state_time = time.time() - res_start

    res_start = time.time()
    optimization_state_impact_scores = attractors. \
        vertex_state_impact_scores(graph_copy, current_attractors=current_attractors,
                                   max_transient_len=optimization_max_transient_len, verbose=False,
                                   relative_attractor_basin_sizes=basin_sizes,
                                   key_slice_size=15)
    optimization_state_time = time.time() - res_start

    result = VertexImpactResult(graph_name=name, is_random=(randomize_functions and randomize_edges),
                                size=graph_name_to_attributes[name]['size'],
                                maximal_change_bits=1,  # TODO: repeat with more?
                                n_inputs=graph_name_to_attributes[name]['n_inputs'],
                                normalized_n_inputs=graph_name_to_attributes[name]['normalized_n_inputs'],
                                max_degree=graph_name_to_attributes[name]['max_degree'],
                                mean_degree=graph_name_to_attributes[name]['mean_degree'],
                                num_attractors=num_attractors,
                                mean_attractor_length=mean_attractor_length,
                                median_attractor_length=median_attractor_length,
                                std_attractor_length=std_attractor_length,
                                max_attractor_length=max_attractor_length,
                                std_attractor_basin=std_attractor_basin,
                                median_attractor_basin=median_attractor_basin,
                                max_attractor_basin=max_attractor_basin,
                                optimization_model_impact_scores=optimization_model_impact_scores,
                                stochastic_model_impact_scores=stochastic_model_impact_scores,
                                optimization_model_addition_impact_scores=optimization_model_addition_impact_scores,
                                stochastic_model_addition_impact_scores=stochastic_model_addition_impact_scores,
                                optimization_state_impact_scores=optimization_state_impact_scores,
                                stochastic_state_impact_scores=stochastic_state_impact_scores,
                                optimization_model_time=optimization_model_time,
                                stochastic_model_time=stochastic_model_time,
                                optimization_model_addition_time=optimization_model_addition_time,
                                stochastic_model_addition_time=stochastic_model_addition_time,
                                optimization_state_time=optimization_state_time,
                                stochastic_state_time=stochastic_state_time
                                )
    print "time taken for graph {} impact scores function: {:.2f} secs".format(name, time.time() - start)
    return result


def main():
    print("stochastic_n_iter={}, parallel={}, n_processes={}, timeout_seconds={}, n_tests={},"
          "filter_out_timed_out_graphs={}, graph_parent_dir={}, optimization_max_len={}, "
          "optimization_max_num={},graph_size_filter={}, queue_all_tasks={},"
          "n_attractors_filter={}".format(stochastic_n_iter, parallel, n_processes,
                                               timeout_seconds, n_tests, filter_out_timed_out_graphs,
                                               graph_parent_dir, optimization_max_len,
                                               optimization_max_num, graph_size_filter, queue_all_tasks,
                                          n_attractors_filter))

    # TODO: use grownups' argument parsing library.
    if len(sys.argv) == 1:
        output_path = 'optimization_and_stochastic_impact_scores.csv'
    else:
        output_path = sys.argv[1]
    print "saving output to {}".format(output_path)

    results = []

    biological_graphs = []
    biological_graph_names = []
    candidate_biological_graph_names = os.listdir(graph_parent_dir)
    for graph_dir in candidate_biological_graph_names:
        try:
            G = graphs.Network.parse_boolean_tables(os.path.join(graph_parent_dir, graph_dir))
            if len(G.vertices) <= graph_size_filter:
                n_estimated_attractors = len(stochastic.estimate_attractors(G, n_walks=stochastic_n_iter, max_walk_len=None,
                                                                        with_basins=False))
                if n_estimated_attractors <= n_attractors_filter:
                    biological_graphs.append(G)
                    biological_graph_names.append(graph_dir)

        except ValueError as e:
            if e.message.startswith("Model export from cellcollective failed"):
                print "warning - did not load graph {}".format(graph_dir)

    print "filtering left {} graphs".format(len(biological_graphs))

    graph_name_to_attributes = dict()
    for i, graph, name in zip(range(len(biological_graphs)), biological_graphs, biological_graph_names):
        n_inputs = len([v for v in graph.vertices if len(v.predecessors()) == 0])
        max_degree = max([len(v.predecessors()) for v in graph.vertices])
        size = len(graph.vertices)
        mean_degree = sum([len(v.predecessors()) for v in graph.vertices]) / float(size)
        normaliezd_n_inputs = n_inputs / float(size)
        graph_name_to_attributes[name] = dict(n_inputs=n_inputs, max_degree=max_degree,
                                                    size=size, mean_degree=mean_degree,
                                              normalized_n_inputs=normaliezd_n_inputs)
        print "#{}; {} input nodes for graph {} of size {} and max degree {}".format(i, n_inputs, name,
                                                                                     size, max_degree)
    results = []
    for test in (range(n_tests) if queue_all_tasks else range(1)):
        test_start = time.time()
        print "test #{}".format(test)
        is_biological = (not only_random) and (test == 0)

        if parallel:
            os.environ['GRB_LICENSE_FILE'] = \
                '/a/home/cc/students/cs/arielbro/{}_gurobi_license/gurobi.lic'.format(platform.node())
            start = time.time()
            repeat = len(biological_graphs)
            # with ProcessPool() as pool:
            pool = multiprocessing.Pool(processes=n_processes)
            future = pool.map(one_graph_impact_score_estimation_wrapper,
                               zip(itertools.cycle(biological_graphs), # a bit more memory efficient.
                                   (biological_graph_names * n_tests) if queue_all_tasks else biological_graph_names,
                                   itertools.repeat(is_biological),
                                   itertools.repeat(graph_name_to_attributes)))
                               # timeout=4000)
            pool.close()
            pool.join()
            # results_iterator = future.result()
            results_iterator = iter(future)  # keeping syntax close to Pebble in case I can get it on nova.

            while True:
                try:
                    result = next(results_iterator)
                    results.append(result)
                except StopIteration:
                    break
                except TimeoutError as e:
                    print "Breaking impact score estimation for timeout after {} seconds".\
                        format(timeout_seconds)
                    results.append(None)
                except attractors.TimeoutError as e:
                    print "Breaking impact score estimation (guorbi timeout)".\
                        format(timeout_seconds)
                    results.append(None)

            print "time_taken for impact scores: {:.2f} secs".format(time.time() - start)

        else:
            for graph_index, graph, name in zip(range(len(biological_graphs)), biological_graphs, biological_graph_names):
                print "working on graph #{}".format(graph_index)
                args = graph, name, is_biological, graph_name_to_attributes
                result = one_graph_impact_score_estimation_wrapper(args)
                results.append(result)

        if filter_out_timed_out_graphs:
            #  On first iteration, filter out graphs for which a timeout occurred.
            original_graph_num = len(biological_graphs)
            good_names = set(result.graph_name for result in results if result is not None)
            good_graph_indices = [i for i, name in enumerate(biological_graph_names) if name in good_names]
            biological_graph_names = [name for name in biological_graph_names if name in good_names]
            biological_graphs = [G for i, G in enumerate(biological_graphs) if i in good_graph_indices]
            assert len(biological_graph_names) == len(biological_graphs)
            print "filtered graphs, {} remaining out of {}".format(len(biological_graphs), original_graph_num)

        # save on each iteration, why not
        with open(output_path, 'wb') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["graph_name", "is_random", "size", "maximal_change_bits",
                            "n_inputs", "normalized_n_inputs", "max_degree", "mean_degree",
                             "num_attractors", "mean_attractor_length", "median_attractor_length",
                             "std_attractor_length", "max_attractor_length", "std_attractor_basin",
                             "median_attractor_basin", "max_attractor_basin",
                             "optimization_model_impact_scores", "stochastic_model_impact_scores",
                             "optimization_model_addition_impact_scores",
                             "stochastic_model_addition_impact_scores",
                             "optimization_state_impact_scores", "stochastic_state_impact_scores",
                             "optimization_model_time", "stochastic_model_time",
                             "optimization_model_addition_time", "stochastic_model_addition_time",
                             "optimization_state_time", "stochastic_state_time"
                             ])
            for impact_result in results:
                if impact_result is None:
                    continue
                writer.writerow(impact_result)
        print "time taken for test #{}: {:.2f} secs".format(test, time.time() - test_start)

if __name__ == "__main__":
    main()
