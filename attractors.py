import numpy as np
import multiprocessing
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
from gurobipy import GRB
import stochastic
import subprocess
from enum import Enum
import pickle
import dill
from pathos.multiprocessing import ProcessingPool as Pool

timeout_seconds = int(0.25 * 60 * 60)  # TODO: refactor as arguments.
dubrova_path = "bns"


class ImpactType(Enum):
    Invalidation = 1
    Addition = 2
    Both = 3


class TimeoutError(Exception):
    pass


def find_num_attractors_multistage(G, use_ilp):
    T = 1
    P = 1
    iteration = 1
    while P <= 2 ** len(G.vertices) and T <= 2 ** len(G.vertices):
        print("iteration {}".format(iteration))
        print("P={}, T={}".format(P, T))
        iteration += 1
        start_time = time.time()

        ATTRACTORS = logic.get_attractorlb_lengthub_formula(G, P, T)
        # print(ATTRACTORS)

        if use_ilp:
            model, _ = ilp.logic_to_ilp(ATTRACTORS)
            model.params.LogToConsole = 0
            model.setObjective(0)
            model.optimize()
            # print("# solutions:", model.SolCount)
            sat = model.SolCount >= 1
        else:
            sat = sympy.satisfiable(ATTRACTORS)
            # sat_iterator = sympy.satisfiable(ATTRACTORS, all_models=True)
            # for model in sat_iterator:
            #     print(model)

        # print("time taken for {} check: {:.2f} seconds".format("ILP" if use_ilp else "SAT", time.time() - start_time))
        # print(sat, '\n')
        if sat:
            print("sat")
            P += 1
        else:
            print("not sat")
            T *= 2
        print('')
    P = P - 1
    print("#attractors:{}".format(P))


def find_num_attractors_onestage(G, max_len=None, max_num=None, use_sat=False, verbose=False,
                                 sampling_bounds=None, use_sampling_for_mip_start=True,
                                 simplify_general_boolean=False,
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
                print("Attractor sampling for non-fixed functions currently unsupported. Skipping.")
                sampling_bounds = None
                break
        if not nonfixed:
            sample_start = time.time()
            # simulate graph to obtain some attractors, then feed them as MIP start to the model.
            attractor_basin_tuples = stochastic.estimate_attractors(G, n_walks=attractor_sampling_num_walks,
                                                                    max_walk_len=attractor_sampling_walk_length)
            num_raw_attractors = len(attractor_basin_tuples)
            attractors = [attractor for (attractor, _) in attractor_basin_tuples if len(attractor) <= T]
            # print("sampled {} suitable attractors (and {} total)".format(len(attractors), num_raw_attractors))
            # print(attractors)
            # print("time taken for attractor sampling: {:.2f} seconds".format(time.time() - sample_start))
            if len(attractors) >= P:
                return P
            if not use_sampling_for_mip_start:
                P = P - len(attractors)

    if use_sat:
        ATTRACTORS, active_logic_vars = logic.get_attractors_formula(G, P, T)
        # for arg in ATTRACTORS.args:
        #     print(arg)

        # print("\n")
        # for active_formula in active_formulas:
        #     print(active_formula)
        # print(ATTRACTORS)

        # print("sat:", sympy.satisfiable(ATTRACTORS & active_logic_vars[0]))
        model, formulas_to_variables = ilp.logic_to_ilp(ATTRACTORS)
        active_ilp_vars = [formulas_to_variables[active_logic_var] for active_logic_var in active_logic_vars]
    else:
        model, a_matrix, v_matrix, state_keys, _ = ilp.attractors_ilp_with_keys(G, T, P, simplify_general_boolean=simplify_general_boolean,
                                                                        slice_size=key_slice_size)
        active_attractors = [a_matrix[p, -1] for p in range(a_matrix.shape[0])]
        if sampling_bounds is not None:
            if use_sampling_for_mip_start:
                ilp.set_mip_start(model, v_matrix, active_attractors, attractors)
            else:
                ilp.add_uniqueness_constraints_from_sampled_attractors(model, state_keys, attractors, 2**key_slice_size,
                                                                       active_attractors, "sampling_non_equality")

    model.setObjective(sum(active_attractors), gurobipy.GRB.MAXIMIZE)
    # model.setParam(gurobipy.GRB.Param.NumericFocus, 3)
    # model.setParam(gurobipy.GRB.Param.OptimalityTol, 1e-6) # gurobi warns against using those for numerical issues
    # model.setParam(gurobipy.GRB.Param.IntFeasTol, 1e-9)
    # model.setParam(gurobipy.GRB.Param.MIPGapAbs, 0.1)
    if not verbose:
        model.params.LogToConsole = 0

    # model.tune()  # try automatic parameter tuning
    # model.getTuneResult(0)  # take best tuning parameters
    # model.write('tune v-{} P-{} T-{}.prm'.format(len(G.vertices), P, T))
    # print(model)
    # model_vars = model.getVars()
    # ilp.print_model_constraints(model)
    # for var in model.getVars():
    #     var.Start = 0

    model.Params.timelimit = timeout_seconds

    model.optimize()
    # model.update()
    # model.write("./model_mip_start.mst") # that just write sthe final solution as a MIP start...

    # ilp.print_opt_solution(model)
    # print(model)
    if model.Status != gurobipy.GRB.OPTIMAL:
        print("warning, model not solved to optimality.")
        if model.Status == gurobipy.GRB.INFEASIBLE:
            # print("writing IIS data to model_iis.ilp")
            # model.computeIIS()
            # model.write("./model_iis.ilp")
            raise RuntimeError("Gurobi failed to reach optimal solution")
        elif model.Status == gurobipy.GRB.TIME_LIMIT:
            raise TimeoutError("Gurobi failed with time_out")
        else:
            raise ValueError("Gurobi failed, status code - {}".format(model.Status))
    else:
        # print("# attractors = {}".format(model.ObjVal))
        # if model.ObjVal != int(round(model.ObjVal)):
        #     print("warning - model solved with non-integral objective function ({})".format(model.ObjVal))
        # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time))
        # ilp.print_attractors(model)
        # ilp.print_model_values(model)
        # ilp.print_model_constraints(model)
        # model.printStats()
        # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time))
        if (sampling_bounds is None) or use_sampling_for_mip_start:
            return int(round(model.objVal))
        else:
            return int(round(model.objVal)) + len(attractors)

    # for constr in model.getConstrs():
    #     print(constr)
    # print([(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName))
    #        if re.match("a_[0-9]*_[0-9]*", v.varName)]  # abs(v.obj) > 1e-6
    # print([(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName))
    #        if re.match("v_[0-9]*_[0-9]*", v.varName)]


def find_model_bitchange_for_new_attractor(G, max_len, verbose=False, key_slice_size=15, use_dubrova=False,
                                           simplify_boolean=False):
    # TODO: add tests!
    start_time = time.time()

    # first, find original model's attractors.
    for v in G.vertices:
        assert (v.function is not None) or len(v.predecessors()) == 0, "input node or undefined function"
    if use_dubrova:
        attractors = find_attractors_dubrova(G=G, dubrova_path=dubrova_path)
        attractors = [att for att in attractors if len(att) <= max_len]
    else:
        attractors = find_attractors_onestage_enumeration(G=G, max_len=max_len, verbose=verbose,
                                                          simplify_general_boolean=simplify_boolean,
                                                          key_slice_size=key_slice_size)

    model, state_keys, a_matrix, function_bitchange_vars = \
        ilp.bitchange_attractor_ilp_with_keys(G, max_len=max_len, slice_size=key_slice_size)
    ilp.add_model_invariant_uniqueness_constraints(model=model, model_state_keys=state_keys, attractors=attractors,
                                                   upper_bound=2 ** key_slice_size, a_matrix=a_matrix)

    model.addConstr(sum(function_bitchange_vars) >= 1) # So Gurobi doesn't work hard proving this.
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
    # print(model)
    # model_vars = model.getVars()
    # ilp.print_model_constraints(model)
    # for var in model.getVars():
    #     var.Start = 0

    model.Params.timelimit = timeout_seconds

    model.optimize()
    # model.update()
    # model.write("./model_mip_start.mst") # that just write sthe final solution as a MIP start...

    # ilp.print_opt_solution(model)
    # print(model)
    if model.Status != gurobipy.GRB.OPTIMAL:
        print("warning, model not solved to optimality.")
        if model.Status == gurobipy.GRB.INFEASIBLE:
            # print("writing IIS data to model_iis.ilp")
            # model.computeIIS()
            # model.write("./model_iis.ilp")
            return numpy.inf
            # raise RuntimeError("Gurobi failed to reach optimal solution")
        elif model.Status == gurobipy.GRB.TIME_LIMIT:
            raise TimeoutError("Gurobi failed with time_out")
        else:
            raise ValueError("Gurobi failed, status code - {}".format(model.Status))
    else:
        # print("# attractors = {}".format(model.ObjVal))
        if model.ObjVal != int(round(model.ObjVal)):
            print("warning - model solved with non-integral objective function ({})".format(model.ObjVal))
        # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time))
        # ilp.print_attractors(model)
        # ilp.print_model_values(model)
        # ilp.print_model_constraints(model)
        # model.printStats()
        # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time))
        return int(round(model.objVal))

        # for constr in model.getConstrs():
        #     print(constr)
        # print([(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName))
        #        if re.match("a_[0-9]*_[0-9]*", v.varName)]  # abs(v.obj) > 1e-6
        # print([(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName))
        #        if re.match("v_[0-9]*_[0-9]*", v.varName)]


def single_model_bitchange_experiment_wrapper(args):
    return single_model_bitchange_experiment(*args)


def single_model_bitchange_experiment(G, perturbed_lines_dict, current_attractors,
                                      relative_attractor_basin_sizes, impact_type, use_dubrova, cur_dubrova_path,
                                      attractor_estimation_n_iter):
    original_functions = [v.function for v in G.vertices]
    for v_index, perturbed_lines in perturbed_lines_dict.items():
        G.vertices[v_index].function = logic.perturb_line(original_functions[v_index], perturbed_lines)

    score = 0
    if (impact_type == ImpactType.Invalidation) or (impact_type == ImpactType.Both):
        for attractor, basin_size in zip(current_attractors, relative_attractor_basin_sizes):
            is_valid = utility.is_attractor_valid(attractor, G)
            if not is_valid:
                score += basin_size
    if (impact_type == ImpactType.Addition) or (impact_type == ImpactType.Both):
        attractors_start = time.time()
        if use_dubrova:
            new_attractors = find_attractors_dubrova(G, cur_dubrova_path, mutate_input_nodes=True)
        else:
            new_attractors = stochastic.estimate_attractors(G, n_walks=attractor_estimation_n_iter,
                                                            max_walk_len=None,
                                                            with_basins=False)
        # print("time taken to calculate new attractors: {:.2f} secs".format(time.time() - attractors_start))

        # print("current attractors:")
        # print(current_attractors)
        # print("new attractors:")
        # print(new_attractors)

        intersection_size = utility.attractor_lists_intersection_size(current_attractors, new_attractors)
        score += (len(new_attractors) - intersection_size) / float(len(current_attractors))

    if impact_type == ImpactType.Both:
        score /= 2

    for i in range(len(G.vertices)):
        G.vertices[i].function = original_functions[i]

    return score

def stochastic_graph_model_impact_score(G, current_attractors, n_iter=100,
                                        use_dubrova=False,
                                        cur_dubrova_path=dubrova_path,
                                        bits_of_change=1,
                                        relative_attractor_basin_sizes=None,
                                        attractor_estimation_n_iter=1000,
                                        parallel_n_jobs=None,
                                        impact_type=ImpactType.Invalidation):
    """
    Returns the mean impact of for uniformly selected bits_of_change bit changes in the functions of nodes
    (not necessarily same node). If impact_type is "Invalidation", impact is defined as the proportion of graph
    states that will no longer belong to the same attractor in the model, which is equivalent to
    the invalidated attractors weighted by their basin size, which is what is actually computed.
    If no basin sizes are given, attractors are weighted equally.
    If impact_type is "Addition", the new attractors added to the model are counted, equally weighted.
    Note that this can give an impact score greater than 1.
    If impact_type is "both", the average of both modes is returned.
    :param G:
    :param current_attractors:
    :param n_iter:
    :param bits_of_change:
    :param relative_attractor_basin_sizes:
    :return:
    """
    # TODO: consider weighting new attractors by their basin sizes too.
    if relative_attractor_basin_sizes is None:
        relative_attractor_basin_sizes = [1 / float(len(current_attractors))] * len(current_attractors)
    for v in G.vertices:
        assert v.function is not None or len(v.predecessors()) == 0, "input node or undefined function"
    start = time.time()
    original_functions = [v.function for v in G.vertices]

    degrees = [len(v.predecessors()) for v in G.vertices]
    if parallel_n_jobs is not None:
        perturbed_lines_dicts = [utility.choose_k_bits_from_vertex_functions(degrees, bits_of_change) for _ in
                                 range(n_iter)]
        # pool = multiprocessing.Pool(processes=parallel_n_jobs)
        pool = Pool(parallel_n_jobs)

        score_list = pool.map(single_model_bitchange_experiment_wrapper,
                              zip(itertools.repeat(G), perturbed_lines_dicts, itertools.repeat(current_attractors),
                                  itertools.repeat(relative_attractor_basin_sizes), itertools.repeat(impact_type),
                                  itertools.repeat(use_dubrova), itertools.repeat(cur_dubrova_path),
                                  itertools.repeat(attractor_estimation_n_iter)))
        pool.close()
        pool.join()
        score = sum(score_list)
    else:
        score = 0.0
        for iteration in range(n_iter):
            if iteration and not iteration % 25:
                print("iteration #{}".format(iteration))
            iteration_start = time.time()
            perturbed_lines_dict = utility.choose_k_bits_from_vertex_functions(degrees, bits_of_change)

            score += single_model_bitchange_experiment(G, perturbed_lines_dict, current_attractors, relative_attractor_basin_sizes,
                                              impact_type, use_dubrova, cur_dubrova_path, attractor_estimation_n_iter)

    score /= n_iter

    # print("time taken for stochastic model impact scores: {:.2f} secs".format(time.time() - start))
    return score


def stochastic_vertex_model_impact_scores(G, current_attractors, n_iter=100, use_dubrova=False,
                                          cur_dubrova_path=dubrova_path,
                                          bits_of_change=1,
                                          relative_attractor_basin_sizes=None,
                                          attractor_estimation_n_iter=1000,
                                          impact_type=ImpactType.Invalidation,
                                          parallel_n_jobs=None):
    """
    For each vertex in G, returns the mean impact of for uniformly selected bit changes
    in its function. If impact_type is "Invalidation", impact is defined as the proportion of graph
    states that will no longer belong to the same attractor in the model, which is equivalent to
    the invalidated attractors weighted by their basin size, which is what is actually computed.
    If no basin sizes are given, attractors are weighted equally.
    If impact_type is "Addition", the new attractors added to the model are counted, equally weighted.
    Note that this can give an impact score greater than 1.
    If impact_type is "both", the average of both modes is returned.
     :param G:
    :param current_attractors:
    :param n_iter:
    :param bits_of_change:
    :param relative_attractor_basin_sizes:
    :return:
    """
    # TODO: consider weighting new attractors by their basin sizes too.
    if relative_attractor_basin_sizes is None:
        relative_attractor_basin_sizes = [1/float(len(current_attractors))] * len(current_attractors)
    start = time.time()
    for v in G.vertices:
        assert v.function is not None or len(v.predecessors()) == 0, "input node or undefined function"
    vertex_scores = []
    for i, v in enumerate(G.vertices):
        vertex_start = time.time()
        if len(v.predecessors()) == 0:
            vertex_scores.append(np.nan)
        else:
            original_function = v.function

            if parallel_n_jobs is not None:
                truth_table_row_indices = random.sample(list(range(2 ** len(v.predecessors()))), bits_of_change)

                perturbed_lines_dicts = [{i: random.sample(list(range(2 ** len(v.predecessors()))), bits_of_change)}
                                         for _ in range(n_iter)]
                # pool = multiprocessing.Pool(processes=parallel_n_jobs)
                pool = Pool(parallel_n_jobs)

                score_list = pool.map(single_model_bitchange_experiment_wrapper,
                                      zip(itertools.repeat(G), perturbed_lines_dicts,
                                          itertools.repeat(current_attractors),
                                          itertools.repeat(relative_attractor_basin_sizes),
                                          itertools.repeat(impact_type),
                                          itertools.repeat(use_dubrova), itertools.repeat(cur_dubrova_path),
                                          itertools.repeat(attractor_estimation_n_iter)))
                pool.close()
                pool.join()
                score = sum(score_list)
            else:
                score = 0.0
                for iteration in range(n_iter):
                    if iteration and not iteration % 25:
                        print("iteration #{}".format(iteration))

                    truth_table_row_indices = random.sample(list(range(2 ** len(v.predecessors()))), bits_of_change)
                    perturbed_lines_dict = {i: truth_table_row_indices}

                    score += single_model_bitchange_experiment(G, perturbed_lines_dict, current_attractors,
                                                               relative_attractor_basin_sizes,
                                                               impact_type, use_dubrova, cur_dubrova_path,
                                                               attractor_estimation_n_iter)

            # print("time taken for vertex {}: {:.2f} secs. {} out of {}".format(v.name, time.time() - vertex_start, i, len(G.vertices)))

            score /= n_iter
            vertex_scores.append(score)
    # print("time taken for stochastic impact scores: {:.2f} secs".format(time.time() - start))
    return vertex_scores


def find_model_bitchange_probability_for_different_attractors(G, n_iter=100, use_dubrova=True, max_len=15):
    # TODO: implement for other types of functions? Implement less expensively?
    # TODO: write tests
    # TODO: implement individual vertex stability.
    if use_dubrova:
        attractors = find_attractors_dubrova(G, dubrova_path=dubrova_path, mutate_input_nodes=True)
    else:
        attractors = find_attractors_onestage_enumeration(G, max_len=max_len)

    n_changes = 0

    for i in range(n_iter):
        if i and not i % 10:
            print("iteration #{}".format(i))
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
            new_attractors = find_attractors_dubrova(copy, dubrova_path, mutate_input_nodes=True)
        else:
            new_attractors = find_attractors_onestage_enumeration(copy, max_len=15)
        if not utility.attractor_sets_equality(attractors, new_attractors):
            n_changes += 1
    return float(n_changes) / n_iter


def single_state_bitchange_experiment_wrapper(args):
    return single_state_bitchange_experiment(*args)


def single_state_bitchange_experiment(G, state_to_attractor_mapping=None, n_bits=1, bitchange_node_indices=None):
    initial_state = stochastic.random_state(G)
    original_attractor = stochastic.walk_to_attractor(G, initial_state, max_walk=None,
                                                      state_to_attractor_mapping=state_to_attractor_mapping)
    unaltered_state = random.choice(original_attractor)
    perturbed_state = list(unaltered_state)  # so we can perturbe it
    if bitchange_node_indices is not None:
        perturbed_indices = bitchange_node_indices
    else:
        possible_perturbation_indices = [i for i in range(len(G.vertices)) if len(G.vertices[i].predecessors()) > 0]
        perturbed_indices = np.random.choice(possible_perturbation_indices, n_bits, replace=False)
    for index in perturbed_indices:
        perturbed_state[index] = 1 - perturbed_state[index]
    perturbed_state = tuple(perturbed_state)  # so it will be hashable
    perturbed_attractor = stochastic.walk_to_attractor(G, perturbed_state, max_walk=None,
                                                       state_to_attractor_mapping=state_to_attractor_mapping)
    # print "original vs perturbed - "
    # print original_attractor
    # print perturbed_attractor
    return not utility.is_same_attractor(original_attractor, perturbed_attractor)


def stochastic_graph_state_impact_score(G, bits_of_change, parallel_n_jobs=None, n_iter=1000):
    """
    Stochastically estimate the probability for a flip in the state of bits_of_change random bits
    in a network's state to move the network from a certain attractor to another attractor
    (or to the basin of another attractor).
    States to change are sampled by uniform selection then simulation to corresponding attractor.
    Input nodes are not considered as flippable.
    :param G:
    :param n_iter:
    :param parallel:
    :return:
    """
    state_to_attractor_mapping = dict()
    # exploit basin mapping memory
    if parallel_n_jobs is not None:
        # pool = multiprocessing.Pool(processes=parallel_n_jobs)
        pool = Pool(parallel_n_jobs)
        bitchange_results = pool.map(single_state_bitchange_experiment_wrapper,
                                     itertools.repeat((G, None, bits_of_change), n_iter))
        pool.close()
        pool.join()
    else:
        bitchange_results = []
        for _ in range(n_iter):
            bitchange_results.append(single_state_bitchange_experiment(G, state_to_attractor_mapping,
                                                                       n_bits=bits_of_change))
    return sum(bitchange_results) / float(n_iter)


def stochastic_vertex_state_impact_scores(G, n_iter=1000, parallel_n_jobs=None):
    """
    Stochastically estimate the probability for a flip in the state of a each vertex a network to move the network
    from a certain attractor to another attractor (or to the basin of another attractor).
    States to change are sampled by uniform selection then simulation to corresponding attractor.
    Impacts of input nodes are non defined.
    :param G:
    :param n_iter:
    :param parallel:
    :return:
    """
    vertex_scores = []
    state_to_attractor_mapping = dict()
    for vertex_index in range(len(G.vertices)):
        vertex_start = time.time()
        if len(G.vertices[vertex_index].predecessors()) == 0:
            vertex_scores.append(np.nan)
            continue
        # print("working on vertex {} ({} of {})".format(G.vertices[vertex_index].name, vertex_index + 1, len(G.vertices)))

        if parallel_n_jobs is not None:
            # pool = multiprocessing.Pool(processes=parallel_n_jobs)
            pool = Pool(parallel_n_jobs)
            bitchange_results = pool.map(single_state_bitchange_experiment_wrapper,
                                         itertools.repeat((G, None, 1, [vertex_index]), n_iter))
            pool.close()
            pool.join()
        else:
            # exploit basin mapping memory
            bitchange_results = []
            for _ in range(n_iter):
                bitchange_results.append(single_state_bitchange_experiment(G, state_to_attractor_mapping,
                                                                           n_bits=1,
                                                                           bitchange_node_indices=[vertex_index]))
        return sum(bitchange_results) / float(n_iter)

        vertex_scores.append(sum(bitchange_results) / float(n_iter))
        # print("time taken for vertex {}: {:.2f} secs. {} out of {}".format(G.vertices[vertex_index].name, time.time() - vertex_start, vertex_index, len(G.vertices)))

    return vertex_scores


def graph_state_impact_score(G, current_attractors, max_transient_len=30, verbose=True,
                              relative_attractor_basin_sizes=None, key_slice_size=15,
                              maximal_bits_of_change=1, timeout_seconds=None):
    """
    Finds the maximal proportion of attractors, possibly weighted by their basin sizes, that can be
    switched to the basin of another attractor by flipping a set of maximal_bits_of_change vertices.
    The same set is assumed to be flipped throughout attractors, at an arbitrary state within it.
    """
    if len(current_attractors) == 1:
        return 0
    if relative_attractor_basin_sizes is None:
        relative_attractor_basin_sizes = [1 / float(len(current_attractors))] * len(current_attractors)

    model = gurobipy.Model()
    is_input = lambda i: len(G.vertices[i].predecessors()) == 0
    perturbed_vertices_indicators = [model.addVar(vtype=gurobipy.GRB.BINARY, name="perturbed_indicator".format(i)
                                                  ) if not is_input(i) else 0 for i in range(len(G.vertices))]
    model.update()
    model.addConstr(sum(perturbed_vertices_indicators) <= maximal_bits_of_change,
                    name="maximal_bits_of_change_constraint")
    objective = 0
    for attractor_index, relative_basin_size in zip(range(len(current_attractors)), relative_attractor_basin_sizes):
        source_state = [model.addVar(vtype=gurobipy.GRB.BINARY,
                                     name="attractor_{}_first_state_{}".format(attractor_index, i))
                        for i in range(len(G.vertices))]
        perturbed_state = [model.addVar(vtype=gurobipy.GRB.BINARY,
                                     name="attractor_{}_perturbed_state_{}".format(attractor_index, i))
                        for i in range(len(G.vertices))]
        model.update()

        for i, indicator, source_vertex, perturbed_vertex in zip(range(len(G.vertices)),
                                                                 perturbed_vertices_indicators,
                                                                 source_state, perturbed_state):
            # don't allow flips of input nodes
            if is_input(i):
                model.addConstr(source_vertex == perturbed_vertex)
            else:
                # Only need to constrain that a zero indicator enforces (perturbed = source).
                model.addConstr(source_vertex - perturbed_vertex <= indicator, name="attractor_{}_perturbed_indicator_"
                                                                            "constraint_{}>=".format(attractor_index, i))
                model.addConstr(perturbed_vertex - source_vertex <= indicator, name="attractor_{}_perturbed_indicator_"
                                                                           "constraint_<={}".format(attractor_index, i))
        if max_transient_len > 0:
            target_state = [model.addVar(vtype=gurobipy.GRB.BINARY,
                                         name="attractor_{}_second_state_{}".format(attractor_index, i))
                            for i in range(len(G.vertices))]
            model.update()
            ilp.add_path_to_model(G, model, path_len=max_transient_len, first_state_vars=perturbed_state,
                                  last_state_vars=target_state)
        else:
            target_state = perturbed_state
        this_attractor_states = current_attractors[attractor_index]
        other_attractor_states = []
        for other_attractor_index in range(len(current_attractors)):
            if other_attractor_index == attractor_index:
                continue
            other_attractor_states.extend(current_attractors[other_attractor_index])
        first_inclusion_indicator = \
            ilp.add_state_inclusion_indicator(model, first_state=source_state,
                                              second_state_set=this_attractor_states,
                                              slice_size=key_slice_size,
                                              prefix="attractor_{}_first_state".format(attractor_index))
        second_inclusion_indicator = \
            ilp.add_state_inclusion_indicator(model, first_state=target_state,
                                              second_state_set=other_attractor_states,
                                              slice_size=key_slice_size,
                                              prefix="attractor_{}_target_state".format(attractor_index))
        model.addConstr(first_inclusion_indicator == 1, name="first_state_inclusion_constraint")
        objective += relative_basin_size * second_inclusion_indicator
    model.setObjective(objective, sense=GRB.MAXIMIZE)
    # print("second inclusion indicator - {}".format(second_inclusion_indicator))
    if not verbose:
        model.params.LogToConsole = 0
    if timeout_seconds is not None:
        model.Params.timelimit = timeout_seconds
    model.optimize()
    # print("Time taken for ILP solve: {:.2f} (T={}, P={})".format(time.time() - start_time, T, P))
    if model.Status != gurobipy.GRB.OPTIMAL:
        if model.Status == gurobipy.GRB.TIME_LIMIT:
            raise TimeoutError("Gurobi timeout")
        else:
            print("writing IIS data to model_iis.ilp")
            model.computeIIS()
            model.write("model_iis.ilp")
            model.write("model.mps")
            raise ValueError("model not solved to optimality.")
    # print(model.ObjVal)
    # ilp.print_model_values(model)
    # ilp.print_model_constraints(model)
    return model.ObjVal


def vertex_state_impact_scores(G, current_attractors, max_transient_len=30, verbose=True,
                               relative_attractor_basin_sizes=None, key_slice_size=15, timeout_seconds=None):
    """
    For each vertex, finds the proportion of attractors, possibly weighted by their basin sizes, that can be
    switched to the basin of another attractor by flipping that vertex.
    """
    if len(current_attractors) == 1:
        return [0 if len(v.predecessors()) > 0 else np.nan for v in G.vertices]
    if relative_attractor_basin_sizes is None:
        relative_attractor_basin_sizes = [1 / float(len(current_attractors))] * len(current_attractors)
    vertex_scores = []
    for vertex_index in range(len(G.vertices)):
        if len(G.vertices[vertex_index].predecessors()) == 0:
            vertex_scores.append(np.nan)
            continue
        vertex_start = time.time()
        # print("working on vertex {} ({} of {})".format(G.vertices[vertex_index].name, vertex_index + 1, len(G.vertices)))

        score = 0
        for attractor_index in range(len(current_attractors)):
            model = gurobipy.Model()
            source_state = [model.addVar(vtype=gurobipy.GRB.BINARY,
                                         name="attractor_{}_first_state_{}".format(attractor_index, i))
                            for i in range(len(G.vertices))]
            model.update()
            perturbed_state = [1 - var if i == vertex_index else var for i, var in enumerate(source_state)]
            if max_transient_len > 0:
                target_state = [model.addVar(vtype=gurobipy.GRB.BINARY,
                                             name="attractor_{}_second_state_{}".format(attractor_index, i))
                                for i in range(len(G.vertices))]
                model.update()
                ilp.add_path_to_model(G, model, path_len=max_transient_len, first_state_vars=perturbed_state,
                                      last_state_vars=target_state)
            else:
                target_state = perturbed_state
            this_attractor_states = current_attractors[attractor_index]
            other_attractor_states = []
            for other_attractor_index in range(len(current_attractors)):
                if other_attractor_index == attractor_index:
                    continue
                other_attractor_states.extend(current_attractors[other_attractor_index])
            first_inclusion_indicator = \
                ilp.add_state_inclusion_indicator(model, first_state=source_state,
                                                  second_state_set=this_attractor_states,
                                                  slice_size=key_slice_size,
                                                  prefix="attractor_{}_first_state".format(attractor_index))
            second_inclusion_indicator = \
                ilp.add_state_inclusion_indicator(model, first_state=target_state,
                                                  second_state_set=other_attractor_states,
                                                  slice_size=key_slice_size,
                                                  prefix="attractor_{}_target_state".format(attractor_index))
            model.addConstr(first_inclusion_indicator == 1, name="first_state_inclusion_constraint")
            model.setObjective(second_inclusion_indicator, sense=GRB.MAXIMIZE)
            # print("second inclusion indicator - {}".format(second_inclusion_indicator))
            if timeout_seconds is not None:
                model.Params.timelimit = timeout_seconds
            if verbose:
                model.params.LogToConsole = 1
            else:
                model.params.LogToConsole = 0
            model.optimize()
            # print("Time taken for ILP solve: {:.2f} (T={}, P={})".format(time.time() - start_time, T, P))
            if model.Status != gurobipy.GRB.OPTIMAL:
                if model.Status == gurobipy.GRB.TIME_LIMIT:
                    raise TimeoutError("Gurobi timeout")
                else:
                    print("writing IIS data to model_iis.ilp")
                    model.computeIIS()
                    model.write("model_iis.ilp")
                    model.write("model.mps")
                    raise ValueError("model not solved to optimality.")
            elif model.ObjVal != int(round(model.ObjVal)):
                if abs(model.ObjVal - int(round(model.ObjVal))) < 0.0001:
                    print("Warning! Integral model solved with value {}. Ignoring as small numeric error.".format(model.ObjVal))
                else:
                    raise ValueError("model solved with non-integral objective function ({})".format(model.ObjVal))
            else:
                is_destructive = int(round(model.ObjVal))
                print(model.ObjVal)
                ilp.print_model_values(model)
                ilp.print_model_constraints(model)
                score += is_destructive * relative_attractor_basin_sizes[attractor_index]
        vertex_scores.append(score)
        # print("time taken for vertex {}: {:.2f} seconds".format(G.vertices[vertex_index].name, time.time() - vertex_start))
    return vertex_scores


def find_attractors_onestage_enumeration(G, max_len=None, verbose=False, simplify_general_boolean=False,
                                         key_slice_size=15):
    T = 2 ** len(G.vertices) if not max_len else max_len
    start_time = time.time()

    model, activity_matrix, v_matrix, _, _ = ilp.attractors_ilp_with_keys(G, T, 1,
                                                                    simplify_general_boolean=simplify_general_boolean,
                                                                    slice_size=key_slice_size)
    active_attractors = [activity_matrix[p, -1] for p in range(activity_matrix.shape[0])]
    model.addConstr(sum(active_attractors) == 1)
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
    # print(model)
    # model_vars = model.getVars()
    # ilp.print_model_constraints(model)
    # for var in model.getVars():
    #     var.Start = 0

    model.Params.timelimit = timeout_seconds

    model.optimize()
    # model.update()
    # model.write("./model_mip_start.mst") # that just write sthe final solution as a MIP start...

    # ilp.print_opt_solution(model)
    # print(model)
    if model.Status != gurobipy.GRB.OPTIMAL:
        if model.Status == gurobipy.GRB.TIME_LIMIT:
            raise TimeoutError("Gurobi failed with time_out")
        elif model.Status == gurobipy.GRB.INFEASIBLE:
            return []
        else:
            raise RuntimeError("Gurobi failed to reach optimal solution")
    else:
        # print("# attractors = {}".format(model.ObjVal))
        if model.ObjVal != int(round(model.ObjVal)):
            print("warning - model solved with non-integral objective function ({})".format(model.ObjVal))
        # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time))
        # ilp.print_attractors_enumeration(model)
        # ilp.print_model_values(model)
        # ilp.print_model_constraints(model)
        # model.printStats()
        # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time))
        # print("number of attractors: {}".format(model.SolCount))
        return ilp.get_model_attractors(model)

    # for constr in model.getConstrs():
    #     print(constr)
    # print([(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName))
    #        if re.match("a_[0-9]*_[0-9]*", v.varName)]  # abs(v.obj) > 1e-6
    # print([(v.varName, v.X) for v in sorted(model.getVars(), key=lambda var: var.varName))
    #        if re.match("v_[0-9]*_[0-9]*", v.varName)]


def find_min_attractors_model(G, max_len=None, min_attractors=None, key_slice_size=15):
    T = max_len if max_len is not None else 2**len(G.vertices)
    iteration = 1
    p_to_models = dict()
    P = min_attractors if min_attractors is not None else 1
    pool_size = 10  # 2000000000
    while True:
        print("iteration {}".format(iteration))
        print("P={}, T={}".format(P, T))
        iteration += 1
        start_time = time.time()
        model, a_matrix, _, _, _ = ilp.attractors_ilp_with_keys(G, max_len=T, max_num=P, find_full_model=True,
                                                                    model_type_restriction=False,
                                                                    slice_size=key_slice_size)
        active_attractors = [a_matrix[p, -1] for p in range(a_matrix.shape[0])]
        model.params.LogToConsole = 0
        model.setObjective(sum(active_attractors))
        model.addConstr(sum(active_attractors) == P)
        model.params.PoolSolutions = pool_size
        model.params.PoolSearchMode = 2
        model.params.MIPGap = 0
        model.optimize()
        if model.SolCount == pool_size and 2*pool_size < 2000000000:
            print("reached pool size capacity ({}). Doubling capacity".format(pool_size))
            pool_size *= 2
            continue
        elif pool_size == 2*pool_size >= 2000000000:
            print("too much solutions, ignoring this P")
            P += 1
            continue
        function_models = set()
        if model.Status != gurobipy.GRB.OPTIMAL:
            p_to_models[P] = function_models
        else:
            if model.ObjVal != int(round(model.ObjVal)):
                print("warning - model solved with non-integral objective function ({})".format(model.ObjVal))
            # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time))
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
            print("Models with {} attractors: {}".format(P-1, len(selected_models)))
            for boolean_model in selected_models:
                print(boolean_model)
        if model.Status != gurobipy.GRB.OPTIMAL:
            break
        P += 1
        continue
    if P == 2 ** len(G.vertices) + 1 and len(p_to_models[P]) != 0:  # all that remain
        print("Models with {} attractors: {}".format(P, len(function_models)))
        for model in function_models:
            print(model)


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
            model, a_matrix, _, _, _ = ilp.attractors_ilp_with_keys(G, T, P,
                                                                     model_type_restriction=model_type_restriction,
                                                                     simplify_general_boolean=False,
                                                                     slice_size=key_slice_size)
        else:
            model, a_matrix, _ = ilp.direct_graph_to_ilp_classic(G, T, P,
                                                                     model_type_restriction=model_type_restriction,
                                                                     simplify_general_boolean=False,
                                                                    slice_size=key_slice_size)
        active_attractors = [a_matrix[p, -1] for p in range(a_matrix.shape[0])]
        model.setObjective(sum(active_attractors), gurobipy.GRB.MAXIMIZE)
        if not verbose:
            model.params.LogToConsole = 0
        start_time = time.time()
        model.optimize()
        # print("Time taken for ILP solve: {:.2f} (T={}, P={})".format(time.time() - start_time, T, P))
        if model.Status != gurobipy.GRB.OPTIMAL:
            print("warning, model not solved to optimality.")
            # print("writing IIS data to model_iis.ilp")
            # model.computeIIS()
            # model.write("model_iis.ilp")
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
                    # print("found {} attractors so far".format(found_attractors))
                if T >= attractor_length_threshold:
                    break
                T *= 2
    # print("Found maximal model with {} attractors".format(found_attractors))
    function_vars = [var for var in model.getVars() if "f_" in var.VarName
                     or "signs" in var.VarName or "threshold" in var.VarName]
    # print(G)
    # ilp.print_model_values(model, model_vars=function_vars)
    # ilp.print_attractors(model)
    function_vars = [(v.VarName, v.X) for v in function_vars]
    if clean_up:
        del model
    return found_attractors, function_vars


def graph_model_impact_score(G, current_attractors, max_len, max_num,
                              impact_types=ImpactType.Invalidation, verbose=True,
                              relative_attractor_basin_sizes=None,
                              normalize_addition_scores=False,
                              maximal_bits_of_change=1,
                              timeout_seconds=None):
    """
    solves an ILP representing the impact of changing at most maximal_bits_of_change bits in truth tables
    of the model's vertices (not necessarily same ones) on attractors -
    impact is defined as either the addition of new attractors to the model, or rendering present attractors
    invalid, depending on impact_types.
    The score is either the relative number of attractors invalidated/added,
    or if supplied a weighted sum with the invalidated attractors' relative basin sizes
    (given in relative_basin_sizes). If both, divide both by two.
    Ideally, score should be in range [0-1], but since more new attractors can be created than there currently
    are, the score isn't actually bound when impact_types is not invalidation.
    If normalize_addition_scores is True, the addition scores are divided by max_num instead, forcing a score
    in range [0-1], regardless of number of current attractors.
    """
    start = time.time()
    for v in G.vertices:
        assert v.function is not None or len(v.predecessors()) == 0, "input node or undefined function"
    slice_size = 15
    original_functions = [v.function for v in G.vertices]
    for v in G.vertices:
        if len(v.predecessors()) > 0:
            v.function = None

    model, a_matrix, v_matrix, state_keys, vertices_f_vars = ilp.attractors_ilp_with_keys(G,
        max_len=0 if impact_types == ImpactType.Invalidation else max_len,
        max_num=0 if impact_types == ImpactType.Invalidation else max_num)

    for i in range(len(G.vertices)):
        G.vertices[i].function = original_functions[i]

    bits_of_change = 0
    for i in range(len(G.vertices)):
        if(len(G.vertices[i].predecessors()) == 0):
            continue
        truth_table_lines = [G.vertices[i].function(*line_in) for line_in in
                            itertools.product([False, True], repeat=len(G.vertices[i].predecessors()))]
        bits_of_change += sum(1 - f_var if truth_table_line else f_var for
                             (f_var, truth_table_line) in zip(vertices_f_vars[i], truth_table_lines))
    model.addConstr(bits_of_change <= maximal_bits_of_change, name="bits_of_change_constraint")
    objective = 0

    if impact_types == ImpactType.Addition or impact_types == ImpactType.Both:
        ilp.add_model_invariant_uniqueness_constraints(model, state_keys, current_attractors,
                                                       upper_bound=2 ** slice_size,
                                                       a_matrix=a_matrix)
        denominator = max_num if normalize_addition_scores else len(current_attractors)
        objective += sum(a_matrix[p, -1] for p in range(max_num)) / denominator

    if impact_types == ImpactType.Invalidation or impact_types == ImpactType.Both:
        if relative_attractor_basin_sizes is None:
            relative_attractor_basin_sizes = \
                [1 / float(len(current_attractors))] * len(current_attractors)
        assert len(relative_attractor_basin_sizes) == len(current_attractors) and \
               abs(sum(relative_attractor_basin_sizes) - 1) < 1e-6, "number of attractors or relative sizes mismatch"
        for attractor, weight in zip(current_attractors, relative_attractor_basin_sizes):
            objective += weight * ilp.add_indicator_for_attractor_invalidity(
                                    model, G, attractor, vertices_f_vars, "invalidated_attractors")
    if impact_types == ImpactType.Both:
        objective /= 2

    # print("time taken to build model: {:.2f}".format(time.time() - start))
    start = time.time()
    if not verbose:
        model.params.LogToConsole = 0
    model.setObjective(objective, gurobipy.GRB.MAXIMIZE)
    if timeout_seconds is not None:
        model.Params.timelimit = timeout_seconds
    model.optimize()
    print "finished optimizing for model impact score"
    if model.Status != gurobipy.GRB.OPTIMAL:
        print("warning, model not solved to optimality.")
        if model.Status == gurobipy.GRB.INFEASIBLE:
            # print("writing IIS data to model_iis.ilp")
            # model.computeIIS()
            # model.write("./model_iis.ilp")
            raise RuntimeError("Gurobi failed to reach optimal solution")
        elif model.Status == gurobipy.GRB.TIME_LIMIT:
            raise TimeoutError("Gurobi failed with time_out")
        else:
            raise ValueError("Gurobi failed, status code - {}".format(model.Status))
    else:
        if model.ObjVal != int(round(model.ObjVal)):
            print("warning - model solved with non-integral objective function ({})".format(model.ObjVal))
        # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start))
    # ilp.print_attractors(model)
    # ilp.print_model_values(model)
    return model.objVal
    # print("score of vertex {}: {:.2f}".format(v.name, vertex_scores[-1]))


def vertex_model_impact_scores(G, current_attractors, max_len, max_num,
                               impact_types=ImpactType.Invalidation, verbose=True,
                               relative_attractor_basin_sizes=None,
                               normalize_addition_scores=False,
                               maximal_bits_of_change=1,
                               timeout_seconds=None):
    """
    For each vertex in G, solves an ILP representing the impact of changing its function on attractors -
    impact is defined as either the addition of new attractors to the model, or rendering present attractors
    invalid, depending on impact_types.
    The number of bits allowed to change in a node's function is given by maximal_bits_of_change. Note that with
    unlimited bits of change, any vertex can invalidate any attractor.
    The score is either the relative number of attractors invalidated/added,
    or if supplied a weighted sum with the invalidated attractors' relative basin sizes
    (given in relative_basin_sizes). If both, divide both by two.
    Ideally, score should be in range [0-1], but since more new attractors can be created than there currently
    are, the score isn't actually bound when impact_types is not invalidation.
    If normalize_addition_scores is True, the addition scores are divided by max_num instead, forcing a score
    in range [0-1], regardless of number of current attractors.
    """
    start = time.time()
    for v in G.vertices:
        assert v.function is not None or len(v.predecessors()) == 0, "input node or undefined function"
    vertex_scores = []
    for i, v in enumerate(G.vertices):
        if len(v.predecessors()) == 0:
            vertex_scores.append(np.nan)
        else:
            slize_size = 15
            original_function = v.function
            v.function = None
            model, a_matrix, v_matrix, state_keys, vertices_f_vars = ilp.attractors_ilp_with_keys(G,
                max_len=0 if impact_types == ImpactType.Invalidation else max_len,
                max_num=0 if impact_types == ImpactType.Invalidation else max_num)
            v.function = original_function
            truth_table_lines = [v.function(*line_in) for line_in in
                                itertools.product([False, True], repeat=len(v.predecessors()))]
            bits_of_change = sum(1 - f_var if truth_table_line else f_var for
                                 (f_var, truth_table_line) in zip(vertices_f_vars[i], truth_table_lines))
            model.addConstr(bits_of_change <= maximal_bits_of_change, name="bits_of_change_constraint")
            objective = 0

            if impact_types == ImpactType.Addition or impact_types == ImpactType.Both:
                ilp.add_model_invariant_uniqueness_constraints(model, state_keys, current_attractors,
                                                               upper_bound=2 ** slize_size,
                                                               a_matrix=a_matrix)
                denominator = max_num if normalize_addition_scores else len(current_attractors)
                objective += sum(a_matrix[p, -1] for p in range(max_num)) / denominator

            if impact_types == ImpactType.Invalidation or impact_types == ImpactType.Both:
                if relative_attractor_basin_sizes is None:
                    relative_attractor_basin_sizes = \
                        [1 / float(len(current_attractors))] * len(current_attractors)
                assert len(relative_attractor_basin_sizes) == len(current_attractors) and \
                       abs(sum(relative_attractor_basin_sizes) - 1) < 1e-6, \
                    "number of attractors or relative size mismatch"
                for attractor, weight in zip(current_attractors, relative_attractor_basin_sizes):
                    objective += weight * ilp.add_indicator_for_attractor_invalidity(
                                            model, G, attractor, vertices_f_vars, "invalidated_attractors")
            if impact_types == ImpactType.Both:
                objective /= 2

            # print("time taken to build model: {:.2f}".format(time.time() - start))
            start = time.time()
            if verbose:
                model.params.LogToConsole = 1
            else:
                model.params.LogToConsole = 0
            model.setObjective(objective, gurobipy.GRB.MAXIMIZE)
            if timeout_seconds is not None:
                model.Params.timelimit = timeout_seconds
            model.optimize()
            if model.Status != gurobipy.GRB.OPTIMAL:
                print("warning, model not solved to optimality.")
                if model.Status == gurobipy.GRB.INFEASIBLE:
                    # print("writing IIS data to model_iis.ilp")
                    # model.computeIIS()
                    # model.write("./model_iis.ilp")
                    raise RuntimeError("Gurobi failed to reach optimal solution")
                elif model.Status == gurobipy.GRB.TIME_LIMIT:
                    raise TimeoutError("Gurobi failed with time_out")
                else:
                    raise ValueError("Gurobi failed, status code - {}".format(model.Status))
            else:
                if model.ObjVal != int(round(model.ObjVal)):
                    print("warning - model solved with non-integral objective function ({})".format(model.ObjVal))
                # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start))
            # ilp.print_attractors(model)
            # ilp.print_model_values(model)
            vertex_scores.append(model.objVal)
            # print("score of vertex {}: {:.2f}".format(v.name, vertex_scores[-1]))
    return vertex_scores


def vertex_degeneracy_scores(G, current_attractors, relative=False, verbose=True):
    """
    For each vertex in G, solves an ILP representing the amount of degeneracy in its function, defined as the
    maximal amount of bits that can be changed in its function without rendering any of the model's current
    attractors invalid.
    If relative is True, divides each score by the number of truth-table rows.
    """
    start = time.time()
    for v in G.vertices:
        assert v.function is not None or len(v.predecessors()) == 0, "input node or undefined function"
    vertex_scores = []
    for i, v in enumerate(G.vertices):
        if len(v.predecessors()) == 0:
            vertex_scores.append(np.nan)
        else:
            slize_size = 15
            original_function = v.function
            v.function = None
            model, a_matrix, v_matrix, state_keys, vertices_f_vars = \
                ilp.attractors_ilp_with_keys(G, max_len=0, max_num=0)
            v.function = original_function
            truth_table_lines = [v.function(*line_in) for line_in in
                                itertools.product([False, True], repeat=len(v.predecessors()))]
            bits_of_change = sum(1 - f_var if truth_table_line else f_var for
                                 (f_var, truth_table_line) in zip(vertices_f_vars[i], truth_table_lines))
            model.setObjective(bits_of_change, gurobipy.GRB.MAXIMIZE)

            for i, attractor in enumerate(current_attractors):
                invalidity_indicator = ilp.add_indicator_for_attractor_invalidity(
                                        model, G, attractor, vertices_f_vars, "invalidated_attractors")
                model.addConstr(invalidity_indicator == 0, name="validity_constraint_attractor_{}".format(i))

            # print("time taken to build model: {:.2f}".format(time.time() - start))
            start = time.time()
            if not verbose:
                model.params.LogToConsole = 0
            model.Params.timelimit = timeout_seconds
            model.optimize()
            if model.Status != gurobipy.GRB.OPTIMAL:
                print("warning, model not solved to optimality.")
                if model.Status == gurobipy.GRB.INFEASIBLE:
                    # print("writing IIS data to model_iis.ilp")
                    # model.computeIIS()
                    # model.write("./model_iis.ilp")
                    raise RuntimeError("Gurobi failed to reach optimal solution")
                elif model.Status == gurobipy.GRB.TIME_LIMIT:
                    raise TimeoutError("Gurobi failed with time_out")
                else:
                    raise ValueError("Gurobi failed, status code - {}".format(model.Status))
            else:
                if model.ObjVal != int(round(model.ObjVal)):
                    print("warning - model solved with non-integral objective function ({})".format(model.ObjVal))
                # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start))
            # ilp.print_attractors(model)
            # ilp.print_model_values(model)
            denominator = 1 if not relative else 2 ** len(v.predecessors())
            vertex_scores.append(model.objVal / denominator)
            # print("score of vertex {}: {:.2f}".format(v.name, vertex_scores[-1]))
    return vertex_scores


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
    print("Time taken:{:.2f} seconds".format(end - start))
    print("Estimated attractors:{}.\nAverage length:{:.2f},\nAverage Basin length found:{:.2f},\nAverage Basin length normalized by coverage: {:.2f}".format(len(attractors), average_length, average_basin, average_basin / coverage_ratio))


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
        print("done {} graphs".format(i + 1))
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
        print("done {} graphs".format(i + 1))
    with open(path, 'w') as output_file:
        writer = csv.writer(output_file)
        writer.writerows(res)


def find_attractors_dubrova(G, dubrova_path, mutate_input_nodes=False):
    """
    Export G, call dubrova's algorithm on it, and parse the attractors from results.
    If mutate_input_nodes is False, assumes all input nodes without a function are constantly off,
    otherwise iterates over value combinations of all None functioned input nodes.
    :param G:
    :param dubrova_path:
    :return: n_attractors
    """
    # TODO: write tests.
    G = G.copy()
    input_nodes = [v for v in G.vertices if len(v.predecessors()) == 0]
    unspecified_input_nodes = [v for v in input_nodes if v.function is None]
    input_node_orig_functions = [v.function for v in input_nodes]
    attractors = []

    for v in G.vertices:
        if len(v.predecessors()) != 0 and v.function is None:
            raise ValueError("Can't run dubrova with a non-fixed vertex function")

    for input_combination in itertools.product([0, 1],
                                               repeat=len(unspecified_input_nodes) if mutate_input_nodes else 0):
        for node, orig_func in zip(input_nodes, input_node_orig_functions):
            if node not in unspecified_input_nodes:
                node.function = bool(node.function()) if callable(node.function) else bool(node.function)
        if mutate_input_nodes:
            for i, input_node in enumerate(unspecified_input_nodes):
                input_node.function = bool(input_combination[i])

        temp_network_path = "./temp_network_{}_{}.cnet".format(int(time.time()), os.getpid())
        graphs.Network.export_to_cnet(G, temp_network_path)
        env = os.environ.copy()
        env['PATH'] += ";C:/cygwin/bin"  # TODO: less hardcoding (it somehow didn't have the right PATH)
        # print("calling dubrova with args={}".format([dubrova_path, temp_network_path]))
        process = subprocess.Popen(args=[dubrova_path, temp_network_path],
                                   stderr=subprocess.STDOUT, stdout=subprocess.PIPE, env=env)
        out, _ = process.communicate()
        if process.returncode != 0:
            raise RuntimeError("Error while running Dubrova - code={}, message={}".
                               format(process.returncode, out))
        os.remove(temp_network_path)
        # Dubrova's output format has a total of attractors on the one before last line.
        # Attractor lengths are last word of lines starting with "Attractor"
        num_attractors = int(out.split(b"\n")[-2].split(b" ")[-1])

        cur_attractor = []
        for line in out.split(b"\n"):
            if line.startswith(b"0") or line.startswith(b"1"):
                cur_attractor.append([int(c) for c in line if c != b" "])
            elif len(cur_attractor) > 0:
                attractors.append(cur_attractor[::-1])  # Appearently Dubrova writes attractors reversed! (bottom->top)
                cur_attractor = []

    return attractors


def find_num_steady_states(G, verbose=False, simplify_general_boolean=False):
    model, v_vars_dict = ilp.steady_state_ilp(G, simplify_general_boolean)
    if not verbose:
        model.params.LogToConsole = 0
    model.params.PoolSolutions = 2000000000
    model.params.PoolSearchMode = 2
    # print(model)
    # model_vars = model.getVars()
    # ilp.print_model_constraints(model)
    start = time.time()
    model.optimize()
    # model.update()

    # ilp.print_opt_solution(model)
    # print(model)
    if model.Status != gurobipy.GRB.OPTIMAL:
        # print("warning, model not solved to optimality.")
        # print("writing IIS data to model_iis.ilp")
        # model.computeIIS()
        # model.write("./model_iis.ilp")
        # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start))
        return 0
    else:
        # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start_time))
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
        # print("steady states:")
        # for ss in steady_states:
        #     print(reduce(lambda x, y: str(x) + ", " + str(y), ss))
        # print("time taken for ILP solve: {:.2f} seconds".format(time.time() - start))
        return n_steady_states


def learn_model_from_experiment_agreement(G, experiments, relax_experiments, max_attractor_len, timeout_seconds=None,
                                          key_slice_size=15, verbose=True, allow_suboptimal=False):
    """
    Given a set of (partially measured) experiments and a full/partial/unknown model, completes the model
    and the unobserved values in the experiments such that the experiments are as close as possible
    representing states in attractors of the model.
    Args:
        G: A graph with a model composed on it. Model can be partial or entirely missing, represented as
        nodes with a None function field.
        experiments: a set of dictionaries mapping a subset of the indices of vertices in G to their
        measured values. Index is 0 based. Vertices not represented (can vary between experiments) are
        free to take any value.
        relax_experiments: A bool, if true then each bit of the experiment is assigned a variable, and the
        disagreement is the negative of the total number of bits in experiments that are changed subject
        to all experiments being a state in the attractors. If false, then bits of the experiments are
        constants, and the disagreement is the number of experiments that can not be placed in an attractor.
        max_attractor_len: Represents the maximal number of states (including first state) to be considered in an
        attractor. If 1, experiments are taken to be steady-states.
        timeout_seconds: A timeout for the ILP solver, function fails with exception if timeout occurs.
    Returns:
        A copy of the graph a completed model, replacing the None functions of vertices.
        The proportion of agreement with experiments, as a number in range [0,1]
    """
    start_time = time.time()
    n = len(G.vertices)
    G = G.copy()
    model = gurobipy.Model()
    objective = 0
    objective_terms = 0
    is_input = lambda i: len(G.vertices[i].predecessors()) == 0

    # assert input nodes appear in all experiments
    for i in range(n):
        if is_input(i):
            assert all([i in experiment for experiment in experiments])

    # Create experiment states, indexed by [experiment][node]
    if relax_experiments:
        experiment_states = [[model.addVar(vtype=gurobipy.GRB.BINARY, name="experiment_{}_node_{}".format(j, i))
                            for i in range(n)] for j in range(len(experiments))]
        for j, experiment in enumerate(experiments):
            for i, val in experiment.items():
                objective += experiment_states[j][i] if val else (1 - experiment_states[j][i])
                objective_terms += 1
    else:
        experiment_states = [[experiment.get(i, model.addVar(vtype=gurobipy.GRB.BINARY,
                                                           name="experiment_{}_node_{}".format(j, i)))
                            for i in range(n)] for j, experiment in enumerate(experiments)]
    model.update()

    # Create model variables for missing functions.
    v_funcs = [None] * n
    for i in range(n):
        if G.vertices[i].function is not None:
            v_funcs[i] = G.vertices[i].function
        elif len(G.vertices[i].predecessors()) > 0:
            f_vars_dict = dict()
            for var_comb_index, var_combination in enumerate(
                    itertools.product((False, True), repeat=len(G.vertices[i].predecessors()))):
                f_vars_dict[var_combination] = model.addVar(
                    vtype=gurobipy.GRB.BINARY, name="f_{}_{}".format(i, var_comb_index))

            f = lambda *input_combination: f_vars_dict[input_combination]
            v_funcs[i] = f

    model.update()

    # Expect (or require) experiments to belong to attractors
    for exp in range(len(experiments)):
        attractor_states = ilp.add_path_to_model(G, model, max_attractor_len,
                                                 experiment_states[exp],
                                                 last_state_vars=None, v_funcs=v_funcs)
        experiment_validity_indicator = ilp.add_state_inclusion_indicator(
                                            model, experiment_states[exp],
                                            attractor_states, slice_size=key_slice_size,
                                            prefix="experiment_{}_state_inclusion".format(exp),
                                            assume_uniqueness=False)
        if relax_experiments:
            model.addConstr(experiment_validity_indicator == 1, "experiment_{}_validity_enforcement".format(exp))
        else:
            objective += experiment_validity_indicator
            objective_terms += 1

    model.setObjective(objective, sense=GRB.MAXIMIZE)
    if not verbose:
        model.params.LogToConsole = 0
    if timeout_seconds is not None:
        model.Params.timelimit = timeout_seconds
    model.optimize()
    print("Time taken for ILP solve: {:.2f} (T={}, #experiments={})".format(time.time() - start_time,
                                                                            max_attractor_len, len(experiments)))
    if model.Status != gurobipy.GRB.OPTIMAL:
        if model.Status == gurobipy.GRB.TIME_LIMIT and not allow_suboptimal:
            raise TimeoutError("Gurobi timeout")
        else:
            print("writing IIS data to model_iis.ilp")
            model.computeIIS()
            model.write("model_iis.ilp")
            model.write("model.mps")
            raise ValueError("model not solved to optimality.")
    # print(model.ObjVal)
    # ilp.print_model_values(model)
    # ilp.print_model_constraints(model)

    # Update the model found
    for i in range(len(G.vertices)):
        if G.vertices[i].function is None and len(G.vertices[i].predecessors()) > 0:
            outputs = []
            cur_func = v_funcs[i]
            for var_comb_index, var_combination in enumerate(
                    itertools.product((False, True), repeat=len(G.vertices[i].predecessors()))):
                outputs.append(cur_func(*var_combination).X)
            input_names = [v.name for v in G.vertices[i].predecessors()]
            G.vertices[i].function = logic.BooleanSymbolicFunc(input_names=input_names, boolean_outputs=outputs)

    # print("printing constraints")
    # ilp.print_model_constraints(model)
    # print("printing values")
    # ilp.print_opt_solution(model)

    return G, (model.objVal / float(objective_terms)), model.Runtime

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
