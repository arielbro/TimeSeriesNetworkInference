import random
import time
import graphs
import attractors
from shutil import copyfile
from collections import namedtuple
import csv
import stochastic
import numpy
import os
import sys
import stat

if __name__ == "__main__":
    # TODO: use grownups' argument parsing library.
    if len(sys.argv) == 1:
        output_path = 'vertex_impact_results_addition.csv'
    else:
        output_path = sys.argv[1]
    print "saving output to {}".format(output_path)

    # copy Dubrova's executable, to prevent mysterious errors when running multiple processes calling it.
    process_specific_dubrova_path = "temp_{}_{}".format(attractors.dubrova_path, os.getpid())
    copyfile(attractors.dubrova_path, process_specific_dubrova_path)
    os.chmod(process_specific_dubrova_path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
    attractors.dubrova_path = process_specific_dubrova_path

    VertexImpactResult = namedtuple("VertexImpactResult", "graph_name random_functions random_edges size "
                                                          "maximal_change_bits n_inputs normalized_n_inputs "
                                                          "max_degree mean_degree "
                                                          "invalidation_scores "
                                                          # "addition score "
                                                          # "both score "
                                                          )
    results = []

    biological_graphs = []
    biological_graph_names = os.listdir("cellcollective_models")
    for graph_dir in biological_graph_names:
        try:
            G = graphs.Network.parse_boolean_tables(os.path.join("cellcollective_models", graph_dir))
            biological_graphs.append(G)
        except ValueError as e:
            if e.message.startswith("Model export from cellcollective failed"):
                print "warning - did not load graph {}".format(graph_dir)

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

    for test in range(1000):
        test_start = time.time()
        print "test #{}".format(test)
        is_biological = (test == 0)
        for graph_index, graph, name in zip(range(len(biological_graphs)), biological_graphs, biological_graph_names):
            print "working on graph #{}".format(graph_index)
            start = time.time()
            graph_copy = graph.copy()

            randomize_functions, randomize_edges = random.choice([(False, True), (True, False), (True, True)])
            if randomize_functions:
                graph_copy.randomize_functions(preserve_truth_ratio=True)
            if randomize_edges:
                # graph_copy.randomize_incoming_edges()
                start = time.time()
                graph_copy.randomize_edges_by_switching()
                print "time taken for graph randomization={:.2f} secs".format(time.time() - start)
            start = time.time()

            attractor_basin_tuples = stochastic.estimate_attractors(graph, n_walks=1000, max_walk_len=100,
                                                                    with_basins=True)
            current_attractors = [pair[0] for pair in attractor_basin_tuples]
            # TODO: insert number of current attractors, and maybe average length, as attributes to save
            basin_sizes = [pair[1] for pair in attractor_basin_tuples]
            sum_sizes = sum(basin_sizes)
            basin_sizes = [size / float(sum_sizes) for size in basin_sizes]
            print "time taken for attractor estimation={:.2f} secs".format(time.time() - start)
            start = time.time()

            impact_scores = attractors.vertex_impact_scores(graph, current_attractors=current_attractors,
                                                            max_len=12, max_num=5, verbose=False,
                                                            impact_types=attractors.ImpactType.Addition,
                                                            normalize_addition_scores=True,
                                                            relative_attractor_basin_sizes=basin_sizes,
                                                            maximal_bits_of_change=1)

            attributes = graph_name_to_attributes[name]
            result = VertexImpactResult(graph_name=name, random_functions=randomize_functions,
                                        random_edges=randomize_edges,
                                        size=graph_name_to_attributes[name]['size'],
                                        maximal_change_bits=1,
                                        n_inputs=graph_name_to_attributes[name]['n_inputs'],
                                        normalized_n_inputs=graph_name_to_attributes[name]['normalized_n_inputs'],
                                        max_degree=graph_name_to_attributes[name]['max_degree'],
                                        mean_degree=graph_name_to_attributes[name]['mean_degree'],
                                        invalidation_scores=impact_scores)
            results.append(result)

            print "time_taken for impact scores: {:.2f} secs".format(time.time() - start)

        # save on each iteration, why not
        with open(output_path, 'wb') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["graph_name", "random_functions", "random_edges", "size", "maximal_change_bits",
                            "n_inputs", "normalized_n_inputs", "max_degree", "mean_degree",
                            "invalidation_scores",
                                            # "addition score",
                                            # "both score",
                                            ])
            for impact_result in results:
                writer.writerow(impact_result)
        print "time taken for test #{}: {:.2f} secs".format(test, time.time() - test_start)