import itertools
import logic
import graphs
import sympy
import time
import ilp
import random
import csv
import gurobipy
import stochastic


def find_num_attractors_multistage(G, use_ilp):
    T = 1
    P = 1
    iteration = 1
    while P <= 2 ** len(G.vertices) and T <= 2 ** len(G.vertices):
        print "iteration {}".format(iteration)
        print "P={}, T={}".format(P, T)
        iteration += 1
        start_time = time.time()

        ATTRACTORS = logic.get_attractorlb_lengthub_formula(G, P, T)
        # print ATTRACTORS

        if use_ilp:
            model, _ = ilp.logic_to_ilp(ATTRACTORS)
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


def find_num_attractors_onestage(G, max_len=None, max_num=None, use_sat=False, verbose=False):
    T = 2 ** len(G.vertices) if not max_len else max_len
    P = 2 ** len(G.vertices) if not max_num else max_num
    start_time = time.time()

    if use_sat:
        ATTRACTORS, active_logic_vars = logic.get_attractors_formula(G, P, T)
        # for arg in ATTRACTORS.args:
        #     print arg

        # print "\n"
        # for active_formula in active_formulas:
        #     print active_formula
        # print ATTRACTORS

        # print "sat:", sympy.satisfiable(ATTRACTORS & active_logic_vars[0])
        model, formulas_to_variables = ilp.logic_to_ilp(ATTRACTORS)
        active_ilp_vars = [formulas_to_variables[active_logic_var] for active_logic_var in active_logic_vars]
    else:
        model, active_ilp_vars = ilp.direct_graph_to_ilp(G, T, P, find_bool_model=False)
    model.setObjective(sum(active_ilp_vars), gurobipy.GRB.MAXIMIZE)
    if not verbose:
        model.params.LogToConsole = 0

    # model.tune()  # try automatic parameter tuning
    # model.getTuneResult(0)  # take best tuning parameters
    # model.write('tune v-{} P-{} T-{}.prm'.format(len(G.vertices), P, T))
    print model
    # for var in model.getVars():
    #     var.Start = 0
    model.optimize()
    # print model
    if model.Status != gurobipy.GRB.OPTIMAL:
        print "warning, model not solved to optimality."
        print "writing IIS data to model_iis.ilp"
        model.computeIIS
        model.write("./model_iis.ilp")

    else:
        print "# attractors = {}".format(model.ObjVal)
    # for constr in model.getConstrs():
    #     print constr
    # print [(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName)
    #        if re.match("a_[0-9]*_[0-9]*", v.varName)]  # abs(v.obj) > 1e-6
    # print [(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName)
    #        if re.match("v_[0-9]*_[0-9]*", v.varName)]
    print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)


def find_min_attractors_model(G):
    # important - destroys G's current model.
    for i, v in enumerate(G.vertices):
        truth_table_size = 2**len(v.predecessors())
        function_variables = [sympy.symbols("f_{}_{}".format(i, j)) for j in range(truth_table_size)]
        v.function = logic.PreRandomizedBooleanSymbolicFunc(function_variables)

    T = 2**len(G.vertices)
    iteration = 1
    p_to_models = dict()
    for P in range(1, 2 ** len(G.vertices) + 1):
        print "iteration {}".format(iteration)
        print "P={}, T={}".format(P, T)
        iteration += 1
        start_time = time.time()
        ATTRACTORS = logic.get_attractorlb_lengthub_formula(G, P, T)
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
    ATTRACTORS = logic.get_attractorlb_lengthub_formula(G, P, T)
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


def stochastic_attractor_estimation(G, n_walks, max_walk_len=None):
    if not max_walk_len:
        max_walk_len = 2**len(G.vertices)

    start = time.time()
    attractors = stochastic.estimate_attractors(G, n_walks=n_walks, max_walk_len=max_walk_len)
    end = time.time()

    total_states = sum([basin for _, basin in attractors])
    average_length = sum(len(attractor) for attractor, _ in attractors) / float(len(attractors))
    average_basin = total_states / float(len(attractors))
    coverage_ratio = float(total_states) / 2**len(G.vertices)
    print "Time taken:{:.2f} seconds".format(end - start)
    print "Estimated attractors:{}.\nAverage length:{:.2f}, " \
          "\nAverage Basin length found:{:.2f}," \
          "\nAverage Basin length normalized by coverage: {:.2f}".format(len(attractors), average_length,
                                                                         average_basin,
                                                                         average_basin / coverage_ratio)


def write_random_graph_estimations_sampling(n_graphs, vertices_bounds, indegree_bounds,
                                            restrict_symmetric_threshold,
                                            restrict_and_or_gates,  n_walks, max_walk_len, path):
    res = [["vertices", "edges", "input_nodes", "attractors", "states_visited", "average_attractor_length", "average_basin_size"]]
    for i in range(n_graphs):
        n = random.randint(*vertices_bounds)
        G = graphs.Network.generate_random(n_vertices=n,
                                       indegree_bounds=indegree_bounds,
                                       restrict_signed_symmetric_threshold=restrict_symmetric_threshold,
                                       restrict_and_or_gates=restrict_and_or_gates)
        input_nodes = len([v for v in G.vertices if len(v.predecessors()) == 0])  # not counting semantic inputs
        attractors = stochastic.estimate_attractors(G, n_walks=min(n_walks, 2*2**n), max_walk_len=max_walk_len)
        total_states = sum([basin for _, basin in attractors])
        average_length = sum(len(attractor) for attractor, _ in attractors) / float(len(attractors))
        average_basin = total_states / float(len(attractors))
        res.append([n, len(G.edges), input_nodes, len(attractors), total_states, average_length, average_basin])
        print "done {} graphs".format(i + 1)
    with open(path, 'w') as output_file:
        writer = csv.writer(output_file)
        writer.writerows(res)

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
# TODO: try and create the ILP without the SAT translation, a lot of the variables might be uneeded.
# TODO:     e.g. a_0_0 >> a_0_1 can just be the constraint a_0_1 >= a_0_0 if it isn't needed in other expressions
