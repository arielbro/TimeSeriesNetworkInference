import utility
import numpy
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

timeout_seconds = 3 * 60  # TODO: refactor somewhere?
dubrova_dir_path = "C:/Users/ariel/Downloads/Attractors - for Ariel/Attractors - for Ariel/BNS_Dubrova_2011"


class TimeoutError(Exception):
    pass


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


def find_num_attractors_onestage(G, max_len=None, max_num=None, use_sat=False, verbose=False,
                                 sampling_bounds=None, use_sampling_for_mip_start=True, simplify_general_boolean=False,
                                 key_slice_size=15):
    # TODO: refactor to a method that returns actual attractors (like in the enumerate version).
    T = 2 ** len(G.vertices) if not max_len else max_len
    P = 2 ** len(G.vertices) if not max_num else max_num

    start_time = time.time()

    if sampling_bounds is not None:
        attractor_sampling_num_walks = sampling_bounds[0]
        attractor_sampling_walk_length = sampling_bounds[1]
        # TODO: implement for graphs with non-fixed vertex functions. Maybe include multiple re-rolls of functions
        # TODO: to find the most attractors.
        nonfixed = False
        for vertex in G.vertices:
            if vertex.function is None and len(vertex.predecessors()) > 0:
                nonfixed = True
                print "Attractor sampling for non-fixed functions currently unsupported. Skipping."
                sampling_bounds = None
                break
        if not nonfixed:
            sample_start = time.time()
            # simulate graph to obtain some attractors, then feed them as MIP start to the model.
            attractor_basin_tuples = stochastic.estimate_attractors(G, n_walks=attractor_sampling_num_walks,
                                                                    max_walk_len=attractor_sampling_walk_length)
            num_raw_attractors = len(attractor_basin_tuples)
            attractors = [attractor for (attractor, _) in attractor_basin_tuples if len(attractor) <= T]
            print "sampled {} suitable attractors (and {} total)".format(len(attractors), num_raw_attractors)
            # print attractors
            print "time taken for attractor sampling: {:.2f} seconds".format(time.time() - sample_start)
            if len(attractors) >= P:
                return P
            if not use_sampling_for_mip_start:
                P = P - len(attractors)

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
        model, active_ilp_vars, v_matrix, state_keys = ilp.attractors_ilp_with_keys(G, T, P, simplify_general_boolean=simplify_general_boolean,
                                                                        slice_size=key_slice_size)
        if sampling_bounds is not None:
            if use_sampling_for_mip_start:
                ilp.set_mip_start(model, v_matrix, active_ilp_vars, attractors)
            else:
                ilp.add_uniqueness_constraints_from_sampled_attractors(model, state_keys, attractors, 2**key_slice_size,
                                                                       active_ilp_vars, "sampling_non_equality")

    model.setObjective(sum(active_ilp_vars), gurobipy.GRB.MAXIMIZE)
    # model.setParam(gurobipy.GRB.Param.NumericFocus, 3)
    # model.setParam(gurobipy.GRB.Param.OptimalityTol, 1e-6) # gurobi warns against using those for numerical issues
    # model.setParam(gurobipy.GRB.Param.IntFeasTol, 1e-9)
    # model.setParam(gurobipy.GRB.Param.MIPGapAbs, 0.1)
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

    model.setParam('TimeLimit', timeout_seconds)

    model.optimize()
    model.update()
    # model.write("./model_mip_start.mst") # that just write sthe final solution as a MIP start...

    # ilp.print_opt_solution(model)
    # print model
    if model.Status != gurobipy.GRB.OPTIMAL:
        print "warning, model not solved to optimality."
        if model.Status == gurobipy.GRB.INFEASIBLE:
            # print "writing IIS data to model_iis.ilp"
            # model.computeIIS()
            # model.write("./model_iis.ilp")
            raise RuntimeError("Gurobi failed to reach optimal solution")
        elif model.Status == gurobipy.GRB.TIME_LIMIT:
            raise TimeoutError("Gurobi failed with time_out")
        else:
            raise ValueError("Gurobi failed, status code - {}".format(model.Status))
    else:
        # print "# attractors = {}".format(model.ObjVal)
        if model.ObjVal != int(round(model.ObjVal)):
            print "warning - model solved with non-integral objective function ({})".format(model.ObjVal)
        # print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)
        # ilp.print_attractors(model)
        # ilp.print_model_values(model)
        # ilp.print_model_constraints(model)
        # model.printStats()
        print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)
        if (sampling_bounds is None) or use_sampling_for_mip_start:
            return int(round(model.objVal))
        else:
            return int(round(model.objVal)) + len(attractors)

    # for constr in model.getConstrs():
    #     print constr
    # print [(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName)
    #        if re.match("a_[0-9]*_[0-9]*", v.varName)]  # abs(v.obj) > 1e-6
    # print [(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName)
    #        if re.match("v_[0-9]*_[0-9]*", v.varName)]


def find_model_bitchange_for_new_attractor(G, max_len, verbose=False, key_slice_size=15, use_dubrova=False,
                                           simplify_boolean=False):
    # TODO: add tests!
    start_time = time.time()

    # first, find original model's attractors.
    for v in G.vertices:
        assert (v.function is not None) or len(v.predecessors()) == 0
    if use_dubrova:
        attractors = find_attractors_dubrova(G=G, dubrova_dir_path=dubrova_dir_path)
        attractors = [att for att in attractors if len(att) <= max_len]
    else:
        attractors = find_attractors_onestage_enumeration(G=G, max_len=max_len, verbose=verbose,
                                                          simplify_general_boolean=simplify_boolean,
                                                          key_slice_size=key_slice_size)

    model, state_keys, a_matrix, function_bitchange_vars = \
        ilp.bitchange_attractor_ilp_with_keys(G, max_len=max_len, slice_size=key_slice_size)
    ilp.add_model_invariant_uniqueness_constraints(model=model, model_state_keys=state_keys, attractors=attractors,
                                                   upper_bound=2 ** key_slice_size, model_activity_vars=a_matrix)

    model.addConstr(sum(function_bitchange_vars) >= 1) # So Gurobi doesn't work hard proving this.
    model.update()
    model.setObjective(sum(function_bitchange_vars), gurobipy.GRB.MINIMIZE)
    # model.setParam(gurobipy.GRB.Param.NumericFocus, 3)
    # model.setParam(gurobipy.GRB.Param.OptimalityTol, 1e-6) # gurobi warns against using those for numerical issues
    # model.setParam(gurobipy.GRB.Param.IntFeasTol, 1e-9)
    # model.setParam(gurobipy.GRB.Param.MIPGapAbs, 0.1)
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

    model.setParam('TimeLimit', timeout_seconds)

    model.optimize()
    model.update()
    # model.write("./model_mip_start.mst") # that just write sthe final solution as a MIP start...

    # ilp.print_opt_solution(model)
    # print model
    if model.Status != gurobipy.GRB.OPTIMAL:
        print "warning, model not solved to optimality."
        if model.Status == gurobipy.GRB.INFEASIBLE:
            # print "writing IIS data to model_iis.ilp"
            # model.computeIIS()
            # model.write("./model_iis.ilp")
            return numpy.inf
            # raise RuntimeError("Gurobi failed to reach optimal solution")
        elif model.Status == gurobipy.GRB.TIME_LIMIT:
            raise TimeoutError("Gurobi failed with time_out")
        else:
            raise ValueError("Gurobi failed, status code - {}".format(model.Status))
    else:
        # print "# attractors = {}".format(model.ObjVal)
        if model.ObjVal != int(round(model.ObjVal)):
            print "warning - model solved with non-integral objective function ({})".format(model.ObjVal)
        # print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)
        # ilp.print_attractors(model)
        # ilp.print_model_values(model)
        # ilp.print_model_constraints(model)
        # model.printStats()
        print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)
        return int(round(model.objVal))

        # for constr in model.getConstrs():
        #     print constr
        # print [(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName)
        #        if re.match("a_[0-9]*_[0-9]*", v.varName)]  # abs(v.obj) > 1e-6
        # print [(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName)
        #        if re.match("v_[0-9]*_[0-9]*", v.varName)]


def find_model_bitchange_probability_for_different_attractors(G, n_iter=100, use_dubrova=True, max_len=15):
    # TODO: implement for other types of functions? Implement less expensively?
    # TODO: write tests
    # TODO: implement individual vertex stability.
    if use_dubrova:
        attractors = find_attractors_dubrova(G, dubrova_dir_path=dubrova_dir_path, mutate_input_nodes=True)
    else:
        attractors = find_attractors_onestage_enumeration(G, max_len=max_len)

    n_changes = 0

    for i in range(n_iter):
        if i and not i % 10:
            print "iteration #{}".format(i)
        copy = G.copy()
        v = random.choice([vertex for vertex in copy.vertices if len(vertex.predecessors()) != 0])
        truth_table_row = random.choice(list(itertools.product([False, True], repeat=len(v.predecessors()))))
        current_value = v.function(*truth_table_row)

        if isinstance(v.function, logic.BooleanSymbolicFunc):
            # don't mutate existing vertex's function
            v.function = logic.BooleanSymbolicFunc(formula=v.function.formula)
            truth_table_row_expression = sympy.And(*[
                var if val else sympy.Not(val) for val, var in zip(truth_table_row, v.function.input_vars)])
            if current_value:
                v.function.formula = sympy.And(v.function.formula, sympy.Not(truth_table_row_expression))
            else:
                v.function.formula = sympy.Or(v.function.formula, truth_table_row_expression)
            truth_table_row_index = sum(2**i for i, val in enumerate(truth_table_row))
            if v.function.boolean_outputs is not None:
                v.function.boolean_outputs[truth_table_row_index] = \
                    not v.function.boolean_outputs[truth_table_row_index]
        elif isinstance(v.function, logic.SymmetricThresholdFunction):
            raise NotImplementedError("can't do a bitchange for a SymmetricThresholdFunction")
        elif callable(v.function):
            old = v.function
            v.function = lambda *state: not current_value if tuple(
                True if p else False for p in state) == truth_table_row else old(*state)
        else:
            raise NotImplementedError("Can't do a bitchange for a non-function")

        if use_dubrova:
            new_attractors = find_attractors_dubrova(copy, dubrova_dir_path, mutate_input_nodes=True)
        else:
            new_attractors = find_attractors_onestage_enumeration(copy, max_len=15)
        if not utility.attractor_sets_equality(attractors, new_attractors):
            n_changes += 1
    return float(n_changes) / n_iter


def find_state_bitchange_probability_for_different_attractors(G, n_iter=1000):
    """
    Stochastically estimate the probability for a bitchange in a network state to move the network
    from a certain attractor to another attractor (or to the basin of another attractor).
    :param G:
    :param n_iter:
    :return:
    """
    # TODO: add an option to distinguish input nodes from others
    # TODO: implement individual attractor stability.
    attractor_changes = 0
    state_to_attractor_mapping = dict()
    for i in range(n_iter):
        # if i and not (i % 10):
        #     print i
        initial_state = stochastic.random_state(G)
        original_attractor = stochastic.walk_to_attractor(G, initial_state, max_walk=None,
                                                          state_to_attractor_mapping=state_to_attractor_mapping)
        unaltered_state = random.choice(original_attractor)
        perturbed_index = random.randint(0, len(unaltered_state) - 1)
        perturbed_state = unaltered_state[:perturbed_index] + \
            (1 - unaltered_state[perturbed_index], ) + unaltered_state[perturbed_index + 1:]
        perturbed_attractor = stochastic.walk_to_attractor(G, perturbed_state, max_walk=None,
                                                           state_to_attractor_mapping=state_to_attractor_mapping)
        if not utility.is_same_attractor(original_attractor, perturbed_attractor):
            attractor_changes += 1
    return attractor_changes / float(n_iter)


def find_attractors_onestage_enumeration(G, max_len=None, verbose=False, simplify_general_boolean=False,
                                         key_slice_size=15):
    T = 2 ** len(G.vertices) if not max_len else max_len
    start_time = time.time()

    model, active_ilp_vars, v_matrix, _ = ilp.attractors_ilp_with_keys(G, T, 1,
                                                                    simplify_general_boolean=simplify_general_boolean,
                                                                    slice_size=key_slice_size)
    model.addConstr(sum(active_ilp_vars) == 1)
    model.params.PoolSolutions = 2000000000
    model.params.PoolSearchMode = 2
    # model.setParam(gurobipy.GRB.Param.NumericFocus, 3)
    # model.setParam(gurobipy.GRB.Param.OptimalityTol, 1e-6) # gurobi warns against using those for numerical issues
    # model.setParam(gurobipy.GRB.Param.IntFeasTol, 1e-9)
    # model.setParam(gurobipy.GRB.Param.MIPGapAbs, 0.1)
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

    model.setParam('TimeLimit', timeout_seconds)

    model.optimize()
    model.update()
    # model.write("./model_mip_start.mst") # that just write sthe final solution as a MIP start...

    # ilp.print_opt_solution(model)
    # print model
    if model.Status != gurobipy.GRB.OPTIMAL:
        if model.Status == gurobipy.GRB.TIME_LIMIT:
            raise TimeoutError("Gurobi failed with time_out")
        elif model.Status == gurobipy.GRB.INFEASIBLE:
            return []
        else:
            raise RuntimeError("Gurobi failed to reach optimal solution")
    else:
        # print "# attractors = {}".format(model.ObjVal)
        if model.ObjVal != int(round(model.ObjVal)):
            print "warning - model solved with non-integral objective function ({})".format(model.ObjVal)
        # print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)
        # ilp.print_attractors_enumeration(model)
        # ilp.print_model_values(model)
        # ilp.print_model_constraints(model)
        # model.printStats()
        print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)
        print "number of attractors: {}".format(model.SolCount)
        return ilp.get_model_attractors(model)

    # for constr in model.getConstrs():
    #     print constr
    # print [(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName)
    #        if re.match("a_[0-9]*_[0-9]*", v.varName)]  # abs(v.obj) > 1e-6
    # print [(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName)
    #        if re.match("v_[0-9]*_[0-9]*", v.varName)]


def find_min_attractors_model(G, max_len=None, min_attractors=None, key_slice_size=15):
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
        model, activity_variables, _, _ = ilp.attractors_ilp_with_keys(G, max_len=T, max_num=P, find_full_model=True,
                                                                    model_type_restriction=False,
                                                                    slice_size=key_slice_size)
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
                             use_state_keys=True, clean_up=False, key_slice_size=15):
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
            model, active_ilp_vars, _, _ = ilp.attractors_ilp_with_keys(G, T, P,
                                                                     model_type_restriction=model_type_restriction,
                                                                     simplify_general_boolean=False,
                                                                     slice_size=key_slice_size)
        else:
            model, active_ilp_vars, _ = ilp.direct_graph_to_ilp_classic(G, T, P,
                                                                     model_type_restriction=model_type_restriction,
                                                                     simplify_general_boolean=False,
                                                                    slice_size=key_slice_size)
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
            model.write("model_iis.ilp")
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


def find_attractors_dubrova(G, dubrova_dir_path, mutate_input_nodes=False):
    """
    Export G, call dubrova's algorithm on it, and parse the attractors from results.
    :param G:
    :param dubrova_dir_path:
    :return: n_attractors
    """
    # TODO: write tests.
    G = G.copy()
    input_nodes = [v for v in G.vertices if len(v.predecessors()) == 0] if mutate_input_nodes else []
    attractors = []

    for v in G.vertices:
        if len(v.predecessors()) != 0 and v.function is None:
            raise ValueError("Can't run dubrova with a non-fixed vertex function")

    for input_combination in itertools.product([0, 1], repeat=len(input_nodes)):

        for i, input_node in enumerate(input_nodes):
            input_node.function = bool(input_combination[i])

        temp_network_path = "./temp_network_{}_{}.cnet".format(int(time.time()), random.randint(1, 1000))
        graphs.Network.export_to_cnet(G, temp_network_path)
        env = os.environ.copy()
        env['PATH'] += ";C:/cygwin/bin"  # TODO: less hardcoding (it somehow didn't have the right PATH)
        process = subprocess.Popen(args=[os.path.join(dubrova_dir_path, "bns.exe"), temp_network_path],
                                   stderr=subprocess.STDOUT, stdout=subprocess.PIPE, env=env)
        out, _ = process.communicate()
        if process.returncode != 0:
            raise RuntimeError("Error while running Dubrova - code={}, message={}".
                               format(process.returncode, out))
        os.remove(temp_network_path)
        # Dubrova's output format has a total of attractors on the one before last line.
        # Attractor lengths are last word of lines starting with "Attractor"
        num_attractors = int(out.split("\n")[-2].split(" ")[-1])

        cur_attractor = []
        for line in out.split("\n"):
            if line.startswith("0") or line.startswith("1"):
                cur_attractor.append([int(c) for c in line if c != " "])
            elif len(cur_attractor) > 0:
                attractors.append(cur_attractor)
                cur_attractor = []

    return attractors


def find_num_steady_states(G, verbose=False, simplify_general_boolean=False):
    model, v_vars_dict = ilp.steady_state_ilp(G, simplify_general_boolean)
    if not verbose:
        model.params.LogToConsole = 0
    model.params.PoolSolutions = 2000000000
    model.params.PoolSearchMode = 2
    # print model
    # model_vars = model.getVars()
    # ilp.print_model_constraints(model)
    start = time.time()
    model.optimize()
    model.update()

    # ilp.print_opt_solution(model)
    # print model
    if model.Status != gurobipy.GRB.OPTIMAL:
        # print "warning, model not solved to optimality."
        # print "writing IIS data to model_iis.ilp"
        # model.computeIIS()
        # model.write("./model_iis.ilp")
        # print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start)
        return 0
    else:
        # print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time)
        # ilp.print_model_values(model)
        # ilp.print_model_constraints(model)
        # model.printStats()
        n_steady_states = model.SolCount
        steady_states = []
        for sol in range(model.SolCount):
            model.setParam(gurobipy.GRB.Param.SolutionNumber, sol)
            steady_state = []
            for v in G.vertices:
                if isinstance(v_vars_dict[v], int):
                    steady_state.append(v_vars_dict[v])
                else:
                    steady_state.append(int(v_vars_dict[v].Xn))
            steady_states.append(steady_state)
        # print "steady states:"
        # for ss in steady_states:
        #     print reduce(lambda x, y: str(x) + ", " + str(y), ss)
        # print "time taken for ILP solve: {:.2f} seconds".format(time.time() - start)
        return n_steady_states



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
