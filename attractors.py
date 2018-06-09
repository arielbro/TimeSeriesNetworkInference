import os
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
import subprocess


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


def find_num_attractors_onestage(G, max_len=None, max_num=None, use_sat=False, verbose=False, require_result=None):
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
        model, active_ilp_vars = ilp.direct_graph_to_ilp_with_keys(G, T, P)
    model.setObjective(sum(active_ilp_vars), gurobipy.GRB.MAXIMIZE)
    model.setParam(gurobipy.GRB.Param.NumericFocus, 3)
    # model.setParam(gurobipy.GRB.Param.OptimalityTol, 1e-6) # gurobi warns against using those for numerical issues
    # model.setParam(gurobipy.GRB.Param.IntFeasTol, 1e-9)
    # model.setParam(gurobipy.GRB.Param.MIPGapAbs, 0.1)
    # TODO: find out why without it I got non-optimal solution values (that I could manually improve without violations)
    if require_result is not None:
        model.addConstr(sum(active_ilp_vars) == require_result, name="optimality_constraint")
        model.update()
    if not verbose:
        model.params.LogToConsole = 0

    # model.tune()  # try automatic parameter tuning
    # model.getTuneResult(0)  # take best tuning parameters
    # model.write('tune v-{} P-{} T-{}.prm'.format(len(G.vertices), P, T))
    # print model
    # model_vars = model.getVars()
    # ilp.print_model_constraints(model)
    # for var in model.getVars():
    #     var.Start = 0
    model.optimize()
    # ilp.print_opt_solution(model)
    # print model
    if model.Status != gurobipy.GRB.OPTIMAL:
        print "warning, model not solved to optimality."
        print "writing IIS data to model_iis.ilp"
        model.computeIIS()
        model.write("./model_iis.ilp")
    else:
        # print "# attractors = {}".format(model.ObjVal)
        if model.ObjVal != int(round(model.ObjVal)):
            print "warning - model solved with non-integral objective function ({})".format(model.ObjVal)
        # print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)
        ilp.print_attractors(model)
        # ilp.print_model_values(model)
        # ilp.print_model_constraints(model)
        # model.printStats()
        print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)
        return int(round(model.objVal))
        # ilp.print_model_values(model, model_vars=model_vars)
    # for constr in model.getConstrs():
    #     print constr
    # print [(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName)
    #        if re.match("a_[0-9]*_[0-9]*", v.varName)]  # abs(v.obj) > 1e-6
    # print [(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName)
    #        if re.match("v_[0-9]*_[0-9]*", v.varName)]


def find_min_attractors_model(G, max_len=None, min_attractors=None):
    T = max_len if max_len is not None else 2**len(G.vertices)
    iteration = 1
    p_to_models = dict()
    P = min_attractors if min_attractors is not None else 1
    pool_size = 10  # 2000000000
    while True:
        print "iteration {}".format(iteration)
        print "P={}, T={}".format(P, T)
        iteration += 1
        start_time = time.time()
        model, activity_variables = ilp.direct_graph_to_ilp_with_keys(G, max_len=T, max_num=P, find_full_model=True,
                                                                      model_type_restriction=
                                                                      graphs.FunctionTypeRestriction.NONE)
        model.params.LogToConsole = 0
        model.setObjective(sum(activity_variables))
        model.addConstr(sum(activity_variables) == P)
        model.params.PoolSolutions = pool_size
        model.params.PoolSearchMode = 2
        model.params.MIPGap = 0
        model.optimize()
        if model.SolCount == pool_size and 2*pool_size < 2000000000:
            print "reached pool size capacity ({}). Doubling capacity".format(pool_size)
            pool_size *= 2
            continue
        elif pool_size == 2*pool_size >= 2000000000:
            print "too much solutions, ignoring this P"
            P += 1
            continue
        function_models = set()
        if model.Status != gurobipy.GRB.OPTIMAL:
            p_to_models[P] = function_models
        else:
            if model.ObjVal != int(round(model.ObjVal)):
                print "warning - model solved with non-integral objective function ({})".format(model.ObjVal)
            # print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)
            for i in range(model.SolCount):
                model.setParam(gurobipy.GRB.Param.SolutionNumber, i)
                function_model = []
                for var in model.getVars():
                    if "f_" in var.VarName:
                        function_model.append((var.VarName, int(round(var.Xn))))
                function_model = tuple(sorted(function_model, key=lambda tup: str(tup[0])))
                function_models.add(function_model)

        p_to_models[P] = function_models
        if P > 1:
            selected_models = [boolean_model for boolean_model in p_to_models[P-1] if boolean_model not in p_to_models[P]]
            print "Models with {} attractors: {}".format(P-1, len(selected_models))
            for boolean_model in selected_models:
                print boolean_model
        if model.Status != gurobipy.GRB.OPTIMAL:
            break
        P += 1
        continue
    if P == 2 ** len(G.vertices) + 1 and len(p_to_models[P]) != 0:  # all that remain
        print "Models with {} attractors: {}".format(P, len(function_models))
        for model in function_models:
            print model


def find_max_attractor_model(G, verbose=False, model_type_restriction=graphs.FunctionTypeRestriction.NONE,
                             attractor_length_threshold=None, attractor_num_threshold=None,
                             use_state_keys=True, clean_up=False):
    """
    Finds a model with maximum attractors for a given graph and function restrictions.
    Performed a modified binary search on sizes of attractors, and a regular one on number of attractors.
    :param G:
    :param verbose:
    :param model_type_restriction: one of none, symmetric-threshold, and-or-gates.
    :param attractor_length_threshold: maximal length of attractors to consider.
    :param attractor_num_threshold: number of attractors to be satisfied with
    :param use_state_keys: whether to use the state-key version of direct graph to ilp.
    :return:
    """
    if attractor_length_threshold is None:
        attractor_length_threshold = 2 ** len(G.vertices)
    if attractor_num_threshold is None:
        attractor_num_threshold = 2 ** len(G.vertices)
    # TODO: consider binary search for T
    P = 1
    T = 1
    while True:
        if use_state_keys:
            model, active_ilp_vars = ilp.direct_graph_to_ilp_with_keys(G, T, P,
                                                                       model_type_restriction=model_type_restriction)
        else:
            model, active_ilp_vars = ilp.direct_graph_to_ilp_classic(G, T, P,
                                                                     model_type_restriction=model_type_restriction)
        model.setObjective(sum(active_ilp_vars), gurobipy.GRB.MAXIMIZE)
        if not verbose:
            model.params.LogToConsole = 0
        start_time = time.time()
        model.optimize()
        # print "Time taken for ILP solve: {:.2f} (T={}, P={})".format(time.time() - start_time, T, P)
        if model.Status != gurobipy.GRB.OPTIMAL:
            print "warning, model not solved to optimality."
            print "writing IIS data to model_iis.ilp"
            model.computeIIS()
            model.write("./model_iis.ilp")
            return None
        else:
            found_attractors = int(round(model.ObjVal))
            if found_attractors >= attractor_num_threshold:
                break
            if found_attractors == P:
                P *= 2
                continue
            else:
                # if found_attractors != P - 1:
                    # print "found {} attractors so far".format(found_attractors)
                if T >= attractor_length_threshold:
                    break
                T *= 2
    # print "Found maximal model with {} attractors".format(found_attractors)
    function_vars = [var for var in model.getVars() if "f_" in var.VarName
                     or "signs" in var.VarName or "threshold" in var.VarName]
    # print G
    # ilp.print_model_values(model, model_vars=function_vars)
    # ilp.print_attractors(model)
    function_vars = [(v.VarName, v.X) for v in function_vars]
    if clean_up:
        del model
    return found_attractors, function_vars


def vertex_impact_scores(G, model_type_restriction=graphs.FunctionTypeRestriction.NONE,
                         attractor_length_threshold=None, attractor_num_threshold=None):
    """
    Iterate over G's vertices, for each one "erease" the knowledge about its function, and find the maximal
    number of attractors using find_max_attractor_model. Score vertices according to that number, and return the
    scores and models achieving them.
    :param G: A graph, assumed to have a complete boolean model defined.
    :param model_type_restriction:
    :param attractor_length_threshold:
    :param attractor_num_threshold:
    :return: vertex_scores[i] = score(v_i), vertex_models[i] = opt f_i.
    """
    for v in G.vertices:
        assert v.function is not None or len(v.predecessors()) == 0
    vertex_scores = []
    vertex_models = []
    for i, v in enumerate(G.vertices):
        # TODO: handle input nodes as a special case.
        # TODO: consider calculating relative impact (ratio between n_attractors with or without the pertubation)
        last_func = v.function
        v.function = None

        # TODO: decide on that or onestage variant
        score, function_var_values = find_max_attractor_model(G, model_type_restriction=model_type_restriction,
                                 attractor_length_threshold=attractor_length_threshold,
                                 attractor_num_threshold=attractor_num_threshold, clean_up=True)
        vertex_scores.append(score)
        vertex_models.append([pair for pair in function_var_values if pair[0].startswith("f_{}_".format(i))])
        print "cur score={}".format(vertex_scores[-1])
        v.function = last_func
    return vertex_scores, vertex_models

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
                                            function_type_restriction,  n_walks, max_walk_len, path):
    res = [["vertices", "edges", "input_nodes", "attractors", "states_visited", "average_attractor_length",
            "average_basin_size"]]
    for i in range(n_graphs):
        n = random.randint(*vertices_bounds)
        G = graphs.Network.generate_random(n_vertices=n,
                                           indegree_bounds=indegree_bounds,
                                           function_type_restriction=function_type_restriction)
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


def write_random_fixed_graph_estimations_sampling(G, n_iter, function_type_restriction, n_walks, max_walk_len, path):
    res = [["vertices", "edges", "input_nodes", "attractors", "states_visited", "average_attractor_length",
            "average_basin_size"]]
    n = len(G.vertices)
    input_nodes = len([v for v in G.vertices if len(v.predecessors()) == 0])  # not counting semantic inputs
    for i in range(n_iter):
        G.randomize_functions(function_type_restriction=function_type_restriction)
        attractors = stochastic.estimate_attractors(G, n_walks=min(n_walks, 2 * 2 ** n), max_walk_len=max_walk_len)
        total_states = sum([basin for _, basin in attractors])
        average_length = sum(len(attractor) for attractor, _ in attractors) / float(len(attractors))
        average_basin = total_states / float(len(attractors))
        res.append([n, len(G.edges), input_nodes, len(attractors), total_states, average_length, average_basin])
        print "done {} graphs".format(i + 1)
    with open(path, 'w') as output_file:
        writer = csv.writer(output_file)
        writer.writerows(res)


def find_num_attractors_dubrova(G, dubrova_dir_path):
    """
    Export G, call dubrova's algorithm on it, and parse number of attractors from results.
    :param G:
    :param dubrova_dir_path:
    :return: n_attractors
    """
    temp_network_path = "./temp_network.cnet"
    graphs.Network.export_to_cnet(G, temp_network_path)
    try:
        return_code = subprocess.call(args=[os.path.join(dubrova_dir_path, "bns"), temp_network_path])
        if return_code >= 0:
            raise Exception("Got an erroneous return code while calling Dubrova - {}".format(return_code))
        else:
            pass
    except Exception as e:
        raise e
    finally:
        os.remove(temp_network_path)


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