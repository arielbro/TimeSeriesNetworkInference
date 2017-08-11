
import itertools
import logic
import graphs
import sympy
import numpy
import functools
import time
import ilp
import random
import csv


def get_attractorlb_lengthub_formula(G, P, T):
    # transform input nodes to stable nodes
    for v in G.vertices:
        if len(v.predecessors()) == 0:
            G.edges.append((v, v))
            v.function = lambda *args: args[0] # functions need to accept a tuple of arguments
    a_matrix = numpy.matrix([[sympy.symbols("a_{}_{}".format(p, t)) for t in range(T+1)] for p in range(P)])
    v_matrix = numpy.array([[[sympy.symbols("v_{}_{}_{}".format(i, p, t)) for t in range(T+1)] for p in range(P)]
                             for i in range(len(G.vertices))])

    ACTIVITY_SWITCH = lambda p: sympy.And(*[a_matrix[p, t] | sympy.And(*[~v_matrix[i, p, t]
                                            for i in range(len(G.vertices))]) for t in range(T)])
    MONOTONE = lambda p: sympy.And(*[a_matrix[p, t] >> a_matrix[p, t+1] for t in range(T)])
    ACTIVE = lambda p: sympy.Or(*[a_matrix[p, t] for t in range(T)])


    predecessors_vars = lambda i, p, t: [v_matrix[vertex.index, p, t] for vertex in G.vertices[i].predecessors()]
    CONSISTENT = lambda p: sympy.And(*[sympy.And(*[
                                     ~a_matrix[p, t] | (sympy.Equivalent(v_matrix[i, p, t+1],
                                     G.vertices[i].function(*predecessors_vars(i, p, t))))
                                     for i in range(len(G.vertices)) if len(G.vertices[i].predecessors()) > 0])
                                     for t in range(T)])
    EQ = lambda p1, p2, t1, t2: sympy.And(*[sympy.Equivalent(v_matrix[i, p1, t1], v_matrix[i, p2, t2])
                                            for i in range(len(G.vertices))])
    CYCLIC = lambda p: (a_matrix[p, 0] & EQ(p, p, 0, T)) | \
                       (sympy.Or(*[~a_matrix[p, t - 1] & a_matrix[p, t] & EQ(p, p, t, T)
                                                                       for t in range(1, T + 1)]))
    SIMPLE = lambda p: sympy.And(*[(a_matrix[p, t] & a_matrix[p, t-1]) >> ~EQ(p, p, t, T) for t in range(1, T)])

    UNIQUE = lambda p1: sympy.And(*[sympy.And(*[a_matrix[p2, t] >> ~EQ(p1, p2, t, T)
                                                for p2 in range(p1 + 1, P)]) for t in range(T + 1)])


    # TODO: write about this in the paper
    STABLE = lambda p: sympy.And(*[sympy.And(*[a_matrix[p, t] & a_matrix[p, t+1] >>
                                               sympy.Equivalent(v_matrix[i, p, t], v_matrix[i, p, t+1])
                                               for t in range(T)]) for i in range(len(G.vertices)) if
                                               len(G.vertices[i].predecessors()) == 0])

    ATTRACTORS = sympy.And(*[MONOTONE(p) & ACTIVE(p) & CONSISTENT(p) & CYCLIC(p) & SIMPLE(p) &
                             UNIQUE(p) & ACTIVITY_SWITCH(p) & STABLE(p) for p in range(P)])

    return ATTRACTORS, a_matrix, v_matrix


def find_num_attractors(G, use_ilp):
    T = 1
    P = 1
    iteration = 1
    while P <= 2 ** len(G.vertices) and T <= 2 ** len(G.vertices):
        print "iteration {}".format(iteration)
        print "P={}, T={}".format(P, T)
        iteration += 1
        start_time = time.time()

        ATTRACTORS, _, _ = get_attractorlb_lengthub_formula(G, P, T)
        # print ATTRACTORS

        if use_ilp:
            model = ilp.logic_to_ilp(ATTRACTORS)
            model.params.LogToConsole = 0
            model.setObjective(0)
            model.optimize()
            # print "# solutions:", model.SolCount
            sat = model.SolCount >= 1
        else:
            sat = sympy.satisfiable(ATTRACTORS)
            # sat_iterator = sympy.satisfiable(ATTRACTORS, all_models=True)
            # for model in sat_iterator:
            #     print model

        print "time taken for {} check: {:.2f} seconds".format("ILP" if use_ilp else "SAT",
                                                               time.time() - start_time)
        # print sat, '\n'
        if sat:
            print "sat"
            P += 1
        else:
            print "not sat"
            T *= 2
        print ''
    P = P - 1
    print "#attractors:{}".format(P)


def find_min_attractors_model(G):
    # important - destroys G's current model.
    for i, v in enumerate(G.vertices):
        truth_table_size = 2**len(v.predecessors())
        function_variables = [sympy.symbols("f_{}_{}".format(i, j)) for j in range(truth_table_size)]
        v.function = functools.partial(logic.variable_bound_boolean_func, function_variables)

    T = 2**len(G.vertices)
    iteration = 1
    p_to_models = dict()
    for P in range(1, 2 ** len(G.vertices) + 1):
        print "iteration {}".format(iteration)
        print "P={}, T={}".format(P, T)
        iteration += 1
        start_time = time.time()
        ATTRACTORS, _, _ = get_attractorlb_lengthub_formula(G, P, T)
        # print ATTRACTORS
        sat_models = sympy.satisfiable(ATTRACTORS, all_models=True)  # TODO: replace with ILP enumeration
        function_models = set()
        for model in sat_models:
            if model:  # take only values of function variables
                function_model = []
                for key, val in model.items():
                    if key.name.startswith("f_"):
                        function_model.append((key, val))
                function_model = tuple(sorted(function_model, key=lambda tup: str(tup[0])))
                function_models.add(function_model)
            else:
                print "not sat!"
        p_to_models[P] = function_models

        if P > 1:
            selected_models = [model for model in p_to_models[P-1] if model not in p_to_models[P]]
            print "Models with {} attractors: {}".format(P-1, len(selected_models))
            for model in selected_models:
                print model
        if P == 2**len(G.vertices):  # all that remain
            print "Models with {} attractors: {}".format(P, len(function_models))
            for model in function_models:
                print model


def sample_graph_sat_results(num_vertices, edge_ratio, T, P):
    start = time.time()
    vertex_names = ["v_" + str(i) for i in range(num_vertices)]
    edges = []
    for possible_edge in itertools.product(vertex_names, vertex_names):
        if random.random() < edge_ratio:
            edges.append(possible_edge)
    G = graphs.Network(vertex_names, edges)
    for v in G.vertices:
        boolean_outputs = [random.choice([True, False]) for i in range(len(v.predecessors()))]
        boolean_function = logic.pre_randomized_boolean_func(boolean_outputs)
        v.function = boolean_function
    ATTRACTORS, _, _ = get_attractorlb_lengthub_formula(G, P, T)
    sat = sympy.satisfiable(ATTRACTORS)
    return ATTRACTORS, time.time() - start


def write_sat_sampling_analysis_table(n_repeats, num_vertices_bound, output_path):
    lines = []
    for i in range(n_repeats):
        n_vertices = random.randint(1, num_vertices_bound)
        edge_ratio = random.random()
        T = random.randint(1, 2**n_vertices)
        P = random.randint(1, 2**n_vertices)
        formula, time_taken = sample_graph_sat_results(n_vertices, edge_ratio, T, P)
        lines.append([n_vertices, edge_ratio, T, P, logic.formula_length(formula), time_taken])
        print "finished iteration #{}".format(i)
    lines.insert(0, ["n_vertices, edge_ratio, T, P, formula_length, time_taken_secs"])
    with open(output_path, 'w') as output_file:
        writer = csv.writer(output_file)
        writer.writerows(lines)

# G = graphs.Network(vertex_names=["A"], edges=[("A", "A")],
#                    vertex_functions=[sympy.Nand])

G = graphs.Network(vertex_names=["A", "B"], edges=[("A", "B"), ("B", "A")],
                   vertex_functions=[sympy.Nand]*2)

# G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "A"), ("C", "A")],
#                    vertex_functions=[sympy.Nand]*3)
#
# G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("B", "C"), ("C", "A")],
#                    vertex_functions=[sympy.Nand]*3)
#
# G = graphs.Network(vertex_names=["A", "B", "C"], edges=[("A", "B"), ("A", "C"), ("B", "C"), ("B", "A"),
#                                                         ("C", "A"), ("C", "B")],
#                    vertex_functions=[sympy.Nor]*3)

# G = graphs.Network.generate_random(8, edge_ratio=0.1)
# print G
# find_num_attractors(G, use_ilp=True)
find_min_attractors_model(G)
# write_sat_sampling_analysis_table(10, 7, "C:/Ariel/Downloads/graph_sampling.csv")


# TODO: think about asynchronous model?
# TODO: problem size analysis.
# TODO: measure and plot running time versus P, T and |V|
# TODO: think/discuss implementing monotone symmetric functions in ILP without the long boolean logic version.
# TODO: discuss minmaxing
# TODO: consider if the activity_switch constraints are wanted
# TODO: For next meeting with Roded, consider given a set of experiments, finding a model that fits best the
# TODO: data, but not considering the data is from a steady-state solution, but from one of the FEW attractors!
# TODO: read MAPK.
# TODO: Ask Hadas about P, T bounds, try to update dynamically (Aracena 2008 has a mathematical, probably NP hard bound)
# TODO: check what SAT solver dubrova used
# TODO: run with reasonable (?) or upper bound T, and increment P. with each iteration, find all
# TODO: possible models (function variable assignments), stop when set changes.
# TODO: create some order on graph states (e.g. lexicographic on the state boolena string), and require each attractor to begin with the lowest value one, and the p'th attractor to start with one less than the p+1'th, so there would be less redundency.

