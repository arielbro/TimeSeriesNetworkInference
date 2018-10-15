import csv
import os
import sys
import graphs
import attractors
import numpy as np
from multiprocessing import Pool
import stochastic

graph_names = []
graph_list = []
is_used = []
for model_directory in ["cellcollective_models", "unused_cellcollective_models/bad_for_modelchange",
                        "unused_cellcollective_models/bad_for_statechange"]:
    current_names = os.listdir(model_directory)
    for graph_dir_name in current_names:
        try:
            G = graphs.Network.parse_boolean_tables(os.path.join(model_directory, graph_dir_name))
            graph_list.append(G)
            graph_names.append(graph_dir_name)
            is_used.append(True if model_directory == "cellcollective_models" else False)
        except ValueError as e:
            if str(e).startswith("Model export from cellcollective failed"):
                print("warning - did not load graph {}".format(graph_dir_name))

print("Finished importing models.")

graph_name_to_attributes = dict()
for i, graph, name in zip(range(len(graph_list)), graph_list, graph_names):
    n_inputs = len([v for v in graph.vertices if len(v.predecessors()) == 0])
    max_degree = max([len(v.predecessors()) for v in graph.vertices])
    size = len(graph.vertices)
    mean_degree = sum([len(v.predecessors()) for v in graph.vertices]) / float(size)
    normaliezd_n_inputs = n_inputs / float(size)
    graph_name_to_attributes[name] = dict(n_inputs=n_inputs, max_degree=max_degree,
                                                size=size, mean_degree=mean_degree,
                                          normalized_n_inputs=normaliezd_n_inputs)
    print("#{}; {} input nodes for graph {} of size {}, mean degree {:.2f} and max degree {}".format(i, n_inputs, name,
                                                                          size, mean_degree, max_degree))

def single_graph_attractors(args):
    is_used, graph = args
    if not is_used:
        return None
    return stochastic.estimate_attractors(G=graph, max_walk_len=500, n_walks=1000, with_basins=True)

pool = Pool(processes=48)
attractor_basin_pair_lists = pool.map(single_graph_attractors, zip(is_used, graph_list))
pool.join()
pool.close()
for graph_name, graph, attractor_basin_pairs in zip(graph_names, graph_list, attractor_basin_pair_lists):
    print(graph_name)
    if attractor_basin_pairs is None:
        continue
    graph_attractors = [pair[0] for pair in attractor_basin_pairs]
    basins = [pair[1] for pair in attractor_basin_pairs]
    graph_name_to_attributes[graph_name]['num_attractors'] = len(graph_attractors)
    graph_name_to_attributes[graph_name]['mean_attractor_len'] = np.mean([len(a) for a in graph_attractors])
    graph_name_to_attributes[graph_name]['max_attractor_len'] = np.max([len(a) for a in graph_attractors])
    graph_name_to_attributes[graph_name]['min_attractor_len'] = np.min([len(a) for a in graph_attractors])
    graph_name_to_attributes[graph_name]['mean_attractor_relative_basin'] = np.mean([len(b) for b in graph_attractors]) / (2 ** len(graph.vertices))
    graph_name_to_attributes[graph_name]['max_attractor_relative_basin'] = np.max([len(b) for b in graph_attractors]) / (2 ** len(graph.vertices))
    graph_name_to_attributes[graph_name]['min_attractor_relative_basin'] = np.min([len(b) for b in graph_attractors]) / (2 ** len(graph.vertices))


with open('graph_properties.csv', 'wb') as f:
    w = csv.DictWriter(f, graph_name_to_attributes[graph_names[0]].keys())
    w.writeheader()
    for attribute_dict in graph_name_to_attributes.values():
        w.writerow(attribute_dict)
