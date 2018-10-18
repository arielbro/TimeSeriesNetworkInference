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
    graph_name_to_attributes[name] = {"name": name, "n inputs": n_inputs, "max degree": max_degree,
                                                "size": size, "mean degree": mean_degree,
                                          "normalized n inputs": normaliezd_n_inputs}
    print("#{}; {} input nodes for graph {} of size {}, mean degree {:.2f} and max degree {}".format(i, n_inputs, name,
                                                                          size, mean_degree, max_degree))

def single_graph_attractors(args):
    is_used, graph = args
    if not is_used:
        return None
    return stochastic.estimate_attractors(G=graph, max_walk_len=500, n_walks=3000, with_basins=True)

def singla_graph_steps_to_attractor(args):
    is_used, graph = args
    if not is_used:
        return None
    return stochastic.estimate_path_len_to_attractor(G=graph, n_iter=3000)

pool = Pool(processes=48)
attractor_basin_pair_lists = pool.map(single_graph_attractors, zip(is_used, graph_list))
steps_to_attractor_lists = pool.map(singla_graph_steps_to_attractor, zip(is_used, graph_list))
pool.close()
pool.join()
for graph_name, graph, attractor_basin_pairs, steps_to_attractor in zip(graph_names, graph_list, attractor_basin_pair_lists, steps_to_attractor_lists):
    print(graph_name)
    if attractor_basin_pairs is None:
        continue
    graph_attractors = [pair[0] for pair in attractor_basin_pairs]
    basins = [pair[1] for pair in attractor_basin_pairs]
    graph_name_to_attributes[graph_name]['num attractors'] = len(graph_attractors)
    graph_name_to_attributes[graph_name]['mean attractor len'] = np.mean([len(a) for a in graph_attractors])
    graph_name_to_attributes[graph_name]['median attractor len'] = np.median([len(a) for a in graph_attractors])
    graph_name_to_attributes[graph_name]['max attractor len'] = np.max([len(a) for a in graph_attractors])
    graph_name_to_attributes[graph_name]['min attractor len'] = np.min([len(a) for a in graph_attractors])
    sum_basins = sum(b for b in basins)
    graph_name_to_attributes[graph_name]['median attractor_relative basin'] = np.median([b for b in basins]) / float(sum_basins)
    graph_name_to_attributes[graph_name]['max attractor relative basin'] = np.max([b for b in basins]) / float(sum_basins)
    graph_name_to_attributes[graph_name]['min attractor relative basin'] = np.min([b for b in basins]) / float(sum_basins)

    graph_name_to_attributes[graph_name]['mean path len to attractor'] = np.mean([p for p in steps_to_attractor])
    graph_name_to_attributes[graph_name]['median path len to attractor'] = np.median([p for p in steps_to_attractor])
    graph_name_to_attributes[graph_name]['min path len to attractor'] = np.min([p for p in steps_to_attractor])
    graph_name_to_attributes[graph_name]['max path len to attractor'] = np.max([p for p in steps_to_attractor])

with open('graph_properties.csv', 'wb') as f:
    w = csv.DictWriter(f, graph_name_to_attributes[graph_names[0]].keys(),
                       delimiter="&", lineterminator="\\\\\n", restval=" ")
    w.writeheader()
    for attribute_dict in graph_name_to_attributes.values():
        w.writerow(attribute_dict)
