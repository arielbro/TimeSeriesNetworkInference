import itertools
import csv
import attractors
import utility
import os
import graphs

"""
A test for the accuracy of model learning under steady-state and attractor assumptions.
Done on specific models, IL-1 and MAPK, for which there are specific changes to the model listed in the
corresponding papers for removing cycles, and GE files I got from Roded.
"""

time_limit = 3600 * 0.4

names = ["IL_1", "EGFR"]
GE_expressions = {"IL_1": utility.parse_ge_file(os.path.join("ge_files", "il1_exp")),
                  "EGFR": utility.parse_ge_file(os.path.join("ge_files", "egfr_exp"))}
acyclic_models = {"IL_1": graphs.Network.parse_boolean_tables(os.path.join("cyclic_acyclic_models", "IL1_acyclic")),
                  "EGFR": graphs.Network.parse_boolean_tables(os.path.join("cyclic_acyclic_models", "EGFR_acyclic"))}
cyclic_models = {"IL_1": graphs.Network.parse_boolean_tables(os.path.join("cyclic_acyclic_models", "IL1_cyclic")),
                  "EGFR": graphs.Network.parse_boolean_tables(os.path.join("cyclic_acyclic_models", "EGFR_cyclic"))}

results = []
parameter_grid = itertools.product(names, [False, True], [False, True], [False, True], [1, 5])
for name, relax_experiments, is_cyclic, is_learnt, T in parameter_grid:
    if not is_cyclic and T != 1:
        continue
    print "name: {}, is_cyclic: {}, is_learnt: {}, relax_experiments: {}, T: {}". \
        format(name, is_cyclic, is_learnt, relax_experiments, T)
    G = cyclic_models[name] if is_cyclic else acyclic_models[name]
    expression_data = GE_expressions[name]
    expression_data = [{G.get_vertex(v_name).index: state for v_name, state in exp.items()} for exp in expression_data]
    if is_learnt:
        G = G.copy()
        for v in G.vertices:
            v.function = None

    _, agreement, runtime = \
        attractors.learn_model_from_experiment_agreement(G, expression_data,
                                                         relax_experiments=relax_experiments,
                                                         max_attractor_len=T,
                                                         timeout_seconds=time_limit,
                                                         allow_suboptimal=True)
    print "name: {}, is_cyclic: {}, is_learnt: {}, relax_experiments: {}, T: {}, agreement: {:.2f}, time: {:.2f}".\
        format(name, is_cyclic, is_learnt, relax_experiments, T, agreement, runtime)
    results.append([name, is_cyclic, is_learnt, T, relax_experiments, agreement, runtime])

    with open("cyclic_learning_results.csv", 'wb') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["model name", "is_cyclic", "is_learnt", "T", "relax_experiments", "agreement", "time"])
        for impact_result in results:
            if impact_result is None:
                continue
            writer.writerow(impact_result)

