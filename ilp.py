import gurobipy
import sympy
import numpy
import itertools
import logic
import time
import math
from graphs import FunctionTypeRestriction
from gurobipy import GRB
from utility import order_key_func

# TODO: find good upper bound again, why didn't 29 work on MAPK_large2?
# http://files.gurobi.com/Numerics.pdf a good resource on numerical issues, high values cause them.


def recursive_logic_to_var(formula, model, formulas_to_variables):
    if formula in formulas_to_variables:
        return formulas_to_variables[formula]
    assert not formula.is_Atom
    arg_vars = [recursive_logic_to_var(arg, model, formulas_to_variables) for arg in formula.args]
    if formula.func == sympy.And:
        new_var = model.addVar(vtype=gurobipy.GRB.BINARY)
        model.update()
        formulas_to_variables[formula] = new_var
        model.addGenConstrAnd(new_var, arg_vars)
        model.update()
        return new_var
    elif formula.func == sympy.Or:
        new_var = model.addVar(vtype=gurobipy.GRB.BINARY)
        formulas_to_variables[formula] = new_var
        model.update()
        model.addGenConstrOr(new_var, arg_vars)
        model.update()
        return new_var
    elif formula.func == sympy.Not:
        # some other constraints need an actual var.
        new_var = model.addVar(vtype=gurobipy.GRB.BINARY)
        formulas_to_variables[formula] = new_var
        model.update()
        model.addConstr(new_var == 1 - arg_vars[0])
        model.update()
        return new_var
    elif formula.func == sympy.Implies:
        new_var = model.addVar(vtype=gurobipy.GRB.BINARY)
        formulas_to_variables[formula] = new_var
        model.update()
        model.addConstr((new_var == 1) >> (arg_vars[1] >= arg_vars[0]))
        model.addConstr((new_var == 0) >> (arg_vars[1] == arg_vars[0] - 1))
        model.update()
        return new_var
    elif formula.func == sympy.Equivalent:
        return recursive_logic_to_var((formula.args[0] >> formula.args[1]) & (formula.args[1] >> formula.args[0]),
                                      model, formulas_to_variables)
        # new_var = model.addVar(vtype=gurobipy.GRB.BINARY)
        # formulas_to_variables[formula] = new_var
        # model.update()
        # model.addConstr((new_var == 1) >> (arg_vars[1] == arg_vars[0]))
        # model.addConstr((new_var == 0) >> (arg_vars[1] != arg_vars[0]))
        # model.update()
        # return new_var
    else:
        raise AssertionError("unexpected formula func encountered: {}".format(formula.func))


def logic_to_ilp(formula):
    symbols = formula.free_symbols
    model = gurobipy.Model()
    formulas_to_variables = {symbol: model.addVar(vtype=gurobipy.GRB.BINARY, name=symbol.name)
                             for symbol in symbols}
    model.update()
    top_formula_var = recursive_logic_to_var(formula, model, formulas_to_variables)  # create ILP vars/constrs
    model.addConstr(top_formula_var == 1, name="top level formula var")  # require top formula to be true
    return model, formulas_to_variables


def unique_state_keys(ordered_state_variables, slice_size):
    """
    Assign a unique numbering to a state, by summing 2**i over active vertices.
    Split the value among several variables, if needed, to fit in 32bit ints with room for mathematical operations.
    :param ordered_state_variables: an ordered fixed iterable of vertex state variables.
    :return: A key, a tuple of integers, identifying this state uniquely among all other states.
    """
    # TODO: see if possible to use larger slices.
    # according to gurobi documentation (https://www.gurobi.com/documentation/7.5/refman/variables.html#subsubsection:IntVars),
    # int values are restrained to +-20 billion, or 2**30.9. In practice, choosing anything close (e.g. 2**28)
    # still results in errors.
    n_parts = int(math.floor(len(ordered_state_variables)/float(slice_size)))
    residue = len(ordered_state_variables) % slice_size
    parts = []
    for i in range(n_parts + 1):
        sum_expression = sum(ordered_state_variables[slice_size*i + j] * 2**j
                             for j in range(slice_size if i != n_parts else residue))
        if i < n_parts or residue != 0:
            parts.append(sum_expression)

    # print "parts={}".format(parts)
    return parts


def create_state_keys_comparison_var(model, first_state_keys, second_state_keys, include_equality, upper_bound,
                                     name_prefix=""):
    """
    Registers and returns a new binary variable, equaling 1 iff first_state_keys > second_state_keys,
    where the order is a lexicographic order (both inputs are tuples of ints, MSB to LSB).
    NOTE: creates len(first_state_keys) - 1 auxiliary variables. Assumes numbers are bounded by upper_bound
    :param first_state_keys:
    :param second_state_keys:
    :param include_equality: boolean, if true then the indicator gets 1 on equality of the two key sets.
    :param name_prefix: a string to prepend to all created variables and constraints
    :return: an indicator variable
    """
    # up to i = n-i
    # Z_i >= (ai - bi + Z_{i + 1}) / (M + 1)
    # Z_i <= (ai - bi + Z_{i + 1} - 1) / (M + 1) + 1
    # Z_n >= (a_n - b_n) / M (can divide by M + 1 for convenience)
    # Z_n <= (a_n - b_n - 1) / (M + 1) + 1
    # multiply by denominator to avoid float inaccuracies
    assert len(first_state_keys) == len(second_state_keys)
    last_var = 0 if not include_equality else 1
    M = upper_bound # M - 1 is actual largest value
    for i in range(len(first_state_keys)):
        a = first_state_keys[-i - 1]
        b = second_state_keys[-i - 1]
        z = model.addVar(vtype=gurobipy.GRB.BINARY, name="{}_{}_indicator".format(name_prefix, i))
        model.update()
        model.addConstr(M * z >= a - b + last_var, name="{}_{}_>=constraint".format(name_prefix, i))
        model.addConstr(M * z <= a - b + last_var + (M - 1), name="{}_{}_<=constraint".format(name_prefix, i))
        last_var = z
        # print "a_{}={}, b_{}={}, M={}, M+1={}".format(len(first_state_keys) - i-1, a, len(first_state_keys) -i-1, b, M, M+1)
    return last_var


def add_uniqueness_constraints_from_sampled_attractors(model, model_state_keys, sampled_attractors, upper_bound,
                                                       last_states_activity_vars, name_prefix):
    """
        Given a set of attractors (ordered sequence of ordered vertex state binary numbers), create constraints
        forcing attractors in a solution to be different than th ones supplied.
        Constraints are made based on the lexicographical order used in the model, so only one state in each sampled
        attractor is compared to one state in each model attractor, and state keys slice size is used for comparisons.
        NOTE: creates (slice_size - 1) * len(sampled_attractors) * P auxiliary variables.
    """
    P = len(last_states_activity_vars)
    T = len(model_state_keys[0]) - 1

    # Record the largest state in each sampled attractor
    largest_sampled_states = list()
    for attractor in sampled_attractors:
        largest_sampled_states.append(max(attractor, key=lambda state: order_key_func(state)))

    # go over pairs of model attractors and sampled attractors and require non-equality.
    slice_size = int(math.log(upper_bound, 2))
    for p, sampled_p in itertools.product(list(range(P)), list(range(len(sampled_attractors)))):
        sampled_state_keys = unique_state_keys(largest_sampled_states[sampled_p], slice_size)
        larger_var = create_state_keys_comparison_var(model, sampled_state_keys, model_state_keys[p][T - 1],
                                                      include_equality=False,
                                                      upper_bound=upper_bound,
                                                      name_prefix="{}_p_{}_sampled_{}_>".format(
                                                          name_prefix, p, sampled_p))
        smaller_var = create_state_keys_comparison_var(model, model_state_keys[p][T - 1], sampled_state_keys,
                                                      include_equality=False,
                                                      upper_bound=upper_bound,
                                                      name_prefix="{}_p_{}_sampled_{}_<".format(
                                                          name_prefix, p, sampled_p))

        # a_p_(T-1) -> larger_var OR smaller_var
        model.addConstr(last_states_activity_vars[p] <= larger_var + smaller_var, name="{}_p_{}_sampled_{}".format(
            name_prefix, p, sampled_p))


def add_indicator_for_attractor_invalidity(model, graph, attractor, vertices_f_vars_list, name_suffix):
    """
        Given an attractor (ordered sequence of ordered vertex state binary numbers) and a list of function truth
         table variables for each node (None for fixed logic nodes), creates an indicator getting a value of
         1 iff the attractor is invalid in the model (i.e. assignment to the function variables).

        NOTE: creates a single indicator variable.
        fixed-logic nodes
    """
    n = len(graph.vertices)
    difference_sum = 0 # accumulate the indicators for node values different than expected in attractor.
    difference_sum_bound = 0
    for step in range(len(attractor)):
        prev_state, cur_state = attractor[step], attractor[(step + 1) % len(attractor)]
        for i in range(n):
            predecessor_indices = tuple(u.index for u in graph.vertices[i].predecessors())
            predecessor_values = tuple(prev_state[index] for index in predecessor_indices)
            if vertices_f_vars_list[i] is None:
                #  It is assumed that the attractor is originally legal, but we'll still assert it
                if len(graph.vertices[i].predecessors()) > 0:
                    assert bool(graph.vertices[i].function(*predecessor_values)) == bool(cur_state[i])
                else:  # input nodes should stay constant
                    assert cur_state[i] == prev_state[i]
            else:
                truth_table_row_index = sum((2**k) * prev_state[prev_index]
                                            for k, prev_index in enumerate(predecessor_indices))
                # create indicator expression for whether model next value is different than attractor value
                model_val = vertices_f_vars_list[i][truth_table_row_index]
                node_difference_indicator = model_val if not cur_state[i] else 1 - model_val
                difference_sum += node_difference_indicator
                difference_sum_bound += 1

    indicator_var = model.addVar(vtype=gurobipy.GRB.BINARY,
                 name="attractor_difference_indicator_{}".format(name_suffix))
    model.update()
    model.addConstr(difference_sum_bound * indicator_var >= difference_sum,
                    name="attractor_difference_indicator_{}_>".format(name_suffix))
    model.addConstr(indicator_var <= difference_sum,
                    name="attractor_difference_indicator_{}_<".format(name_suffix))
    model.update()
    return indicator_var


def add_model_invariant_uniqueness_constraints(
            model, model_state_keys, attractors, upper_bound,
            a_matrix):
    """
        Given a set of attractors (ordered sequence of ordered vertex state binary numbers), create constraints
        forcing attractors in a solution to be different than the ones supplied.
        Constraints are made based on the lexicographical order used in the model, so only one state in each sampled
        attractor is compared to one state in each model attractor, and state keys slice size is used for comparisons.
        NOTE: creates (slice_size - 1) * len(sampled_attractors) * T auxiliary variables.
    """
    if len(attractors) == 0:
        return
    P, T = a_matrix.shape[0], a_matrix.shape[1] - 1  # a_matrix has representation for T + 1 states
    n = len(attractors[0][0])

    def order_key_func(node_states): return sum(node * 2**i for (i, node) in enumerate(node_states))

    # Rotate the supplied attractors to fit the desired order.
    ordered_attractors = []
    for attractor in attractors:
        largest_state_index = max(range(len(attractor)), key=lambda i: order_key_func(attractor[i]))
        ordered_attractors.append([attractor[(t + largest_state_index + 1) % len(attractor)]
                                   for t in range(len(attractor))])
    del attractors
    # go over pairs of model attractors and sampled attractors and require non-equality.
    slice_size = int(math.log(upper_bound, 2))
    for given_p in list(range(len(ordered_attractors))):
        if len(ordered_attractors[given_p]) > T:
            print "warning, given attractor of length {} as reference for a problem with T={}. Ignoring". \
                format(len(ordered_attractors[given_p]), T)
            continue
        for p in list(range(P)):
            # require that either attractor length or actual states be different from the given attractor.

            # length difference (recall model has one redundant step)
            model_attractor_len = sum(a_matrix[p, t] for t in range(a_matrix.shape[1] - 1))
            attractor_len = len(ordered_attractors[given_p])
            # In principle, we want an indicator z for (model_attractor_len != attractor_len).
            # We'll split z to z_+ and z_-, and require only that z_+ ==> model_attractor_len > attractor_len
            #                                               and z_- ==> model_attractor_len < attractor_len
            # Then if inequality holds, (z- + z+) *can* be required to be 1, otherwise both are zero.
            different_length_positive = model.addVar(vtype=gurobipy.GRB.BINARY,
                                                      name="length_difference_{}_{}_>".format(p, given_p))
            different_length_negative = model.addVar(vtype=gurobipy.GRB.BINARY,
                                                      name="length_difference_{}_{}_<".format(p, given_p))
            model.update()
            # with x = (model_attractor_len - attractor_len) in range [-T, T], require (T+1)Z_+ <= x + T,
            model.addConstr((T + 1) * different_length_positive <= (model_attractor_len - attractor_len) + T,
                            name="length_difference_{}_{}_>")
            model.addConstr((T + 1) * different_length_negative <= (attractor_len - model_attractor_len) + T,
                            name="length_difference_{}_{}_<")

            # Now add indicators for difference in corresponding states in each attractor. We want either those to be
            # different, or the length (it is possible for length to differ but states not, because 0^n encodes a real state
            # as well as a vacant state. However, if lengths are the same, there is only one representation and all states
            # should match.

            state_difference_indicators = []
            for t in range(T):  # no need to require difference for last state in model, since it is redundant.
                is_active = t >= T - len(ordered_attractors[given_p])
                if is_active:
                    given_state_keys = unique_state_keys(
                        ordered_attractors[given_p][t - (T - len(ordered_attractors[given_p]))], slice_size)
                    # print "t={}, given_index={}, given_state_keys={}".format(
                    #     t, t - (T - len(ordered_attractors[given_p])), given_state_keys)
                else:
                    given_state_keys = unique_state_keys([0] * n, slice_size)
                larger_var = create_state_keys_comparison_var(model, given_state_keys, model_state_keys[p][t],
                                                              include_equality=False,
                                                              upper_bound=upper_bound,
                                                              name_prefix="p_{}_given_p_{}_t_{}_>".format(p, given_p, t))
                smaller_var = create_state_keys_comparison_var(model, model_state_keys[p][t], given_state_keys,
                                                               include_equality=False,
                                                               upper_bound=upper_bound,
                                                               name_prefix="p_{}_given_p_{}_t_{}__<".format(p, given_p, t))
                state_difference_indicators.extend([larger_var, smaller_var])

            model.update()
            model.addConstr(different_length_positive + different_length_negative +
                            sum(state_difference_indicators) >= 1,
                            name="p_{}_given_p_{}_difference_constraint".format(p, given_p))


def add_truth_table_consistency_constraints(model, v_func, v_next_state_var, predecessors_cur_vars,
                                           name_prefix, activity_variable=None, find_model_f_vars=None):
    """
    Adds consistency constraints to a model, as in Roded's paper.
    :param G:
    :param model:
    :param find_model:
    :return:
    """
    in_degree = len(predecessors_cur_vars)
    for var_comb_index, var_combination in enumerate(itertools.product((False, True), repeat=in_degree)):
        if find_model_f_vars is not None:
            desired_val = find_model_f_vars[var_comb_index]
        else:
            desired_val = 1 if v_func(*var_combination) else 0  # == because sympy
        # this expression is |in_degree| iff their states agrees with var_combination
        indicator_expression = sum(v if state else 1 - v for (v, state) in
                                   zip(predecessors_cur_vars, var_combination))
        # a[p, t] & (indicator_expression = in_degree) => v[i,p,t+1] = f(var_combination).
        # For x&y => a=b, require a <= b + (2 -x -b), a >= b - (2 -x -y)
        activity_expression = activity_variable - 1 if activity_variable is not None else 0
        model.addConstr(v_next_state_var >= desired_val -
                        (in_degree - indicator_expression - activity_expression),
                        name="{}_>=_{}".format(name_prefix, var_comb_index))
        model.addConstr(v_next_state_var <= desired_val +
                        (in_degree - indicator_expression - activity_expression),
                        name="{}_<=_{}".format(name_prefix, var_comb_index))


def build_logic_function_vars(formula, model, name_prefix, symbols_to_variables_dict):
    """
    Builds an indicator for a boolean function recursively.
    :param boolean_function:
    :param model:
    :param formula_to_variables_dict: maps the formula symbols to model variables
    :return:
    """
    # TODO: support more operation types
    if isinstance(formula, sympy.symbol.Symbol):
        return symbols_to_variables_dict[formula]
    elif isinstance(formula, sympy.And):
        argument_vars = [build_logic_function_vars(argument, model, "{}_And_args_{}".format(name_prefix, i),
                                                   symbols_to_variables_dict)
                                                   for (i, argument) in enumerate(formula.args)]
        andVar = model.addVar(vtype=gurobipy.GRB.BINARY,
                              name="{}_And".format(name_prefix))
        for (i, argument_var) in enumerate(argument_vars):
            model.addConstr(andVar <= argument_var, name="{}_And_con_<=_{}".format(name_prefix, i))
        model.addConstr(andVar >= sum(argument_vars) - (len(argument_vars) - 1),
                        name="{}_And_con_>=_{}".format(name_prefix, i))
        return andVar
    elif isinstance(formula, sympy.Or):
        argument_vars = [build_logic_function_vars(argument, model, "{}_Or_args_{}".format(name_prefix, i),
                                                   symbols_to_variables_dict)
                                                   for (i, argument) in enumerate(formula.args)]
        orVar = model.addVar(vtype=gurobipy.GRB.BINARY,
                              name="{}_Or".format(name_prefix))
        for (i, argument_var) in enumerate(argument_vars):
            model.addConstr(orVar >= argument_var, name="{}_Or_con_>=_{}".format(name_prefix, i))
        model.addConstr(orVar <= sum(argument_vars),
                        name="{}_And_res_<=_{}".format(name_prefix, i))
        return orVar
    elif isinstance(formula, sympy.Not):
        try:
            return 1 - symbols_to_variables_dict[formula.args[0]]
        except Exception as e:
            raise e
    elif formula in [sympy.true, sympy.false]:
        return 1 if formula is sympy.true else 0
    else:
        raise NotImplementedError


def add_simplified_consistency_constraints(model, v_func, v_next_state_var, predecessors_cur_vars,
                                           name_prefix, activity_variable=None):
    """
    Adds consistency constraints to the model by transforming a vertex' logic function to an indicator variable,
    and enforcing its value to be equal to the next state's value.
    :param G:
    :param model:
    :param v_next_state_var:
    :param v_cur_var:
    :param predecessors_cur_vars:
    :return:
    """
    if v_func is None:
        raise AssertionError("Nodes with unset functions shouldn't be handled here")
    if isinstance(v_func, sympy.FunctionClass):  # TODO: convert to SymbolicBooleanFunction?
        # instantiate the function expression
        arg_symbols = [sympy.symbols("x_{}".format(j)) for j in range(len(predecessors_cur_vars))]
        func_expression = v_func(*arg_symbols)
    elif isinstance(v_func, logic.BooleanSymbolicFunc):
        func_expression = v_func.formula
        arg_symbols = v_func.input_vars
    else:
        # try a constant function
        try:
            if v_func(None) in [True, False, sympy.true, sympy.false]:
                func_expression = sympy.true if v_func(None) else sympy.false
                arg_symbols = [None] * len(predecessors_cur_vars)
            else:
                raise ValueError("Unkown type of function - " + str(type(v_func)))
        except TypeError:
            raise ValueError("Unkown type of function (non-constant lambda functions aren't allowed)")
    simplified_formula = sympy.simplify(func_expression)
    symbols_to_variables_dict = dict(zip(arg_symbols, predecessors_cur_vars))
    func_var = build_logic_function_vars(simplified_formula, model, name_prefix,
                                         symbols_to_variables_dict)
    if activity_variable is None:
        model.addConstr(v_next_state_var == func_var)
    else:
        activity_part = activity_variable - 1 if activity_variable is not None else 0
        model.addConstr(v_next_state_var >= func_var + activity_part,
                        name="{}_>=".format(name_prefix))
        model.addConstr(v_next_state_var <= func_var - activity_part,
                        name="{}_<=".format(activity_part))


def add_state_inclusion_indicator(model, first_state, second_state_set, slice_size, prefix=None):
    """
    Adds and returns a binary indicator for whether one network state is included in a set of others.
    States should be either iterables of constant binary values (False/True/0/1), or model variables.
    It is assumed that all states in second_state_set are unique
    Indicator is created with use of state hashing and order indicators for state pairs.
    :param model:
    :param first_state:
    :param second_state_set:
    :return:
    """
    indicator_sum = 0
    first_state_keys = unique_state_keys(first_state, slice_size=slice_size)
    # print "len of second state set - {}".format(len(second_state_set))
    for i, second_state in enumerate(second_state_set):
        second_state_keys = unique_state_keys(second_state, slice_size=slice_size)
        larger_var = create_state_keys_comparison_var(model, first_state_keys, second_state_keys,
                                                      include_equality=True,
                                                      upper_bound=2**slice_size,
                                                      name_prefix="{}_state_inclusion_{}_>=".format(prefix, i))
        smaller_var = create_state_keys_comparison_var(model, second_state_keys, first_state_keys,
                                                       include_equality=True,
                                                       upper_bound=2**slice_size,
                                                       name_prefix="{}_state_inclusion_{}_<=".format(prefix, i))
        # print "indicator sum pre: {}".format(indicator_sum)
        indicator_sum += larger_var + smaller_var
        # print "indicator sum post: {}".format(indicator_sum)
    model.update()
    # We want an indicator for inclusion of first_state in the second state, which is equivalent to its equality
    # with exactly one state there (since they're unique), or len(second_state_set) + 1 indicators with value 1.
    # print "indicator sum - {}".format(indicator_sum)
    inclusion_indicator = indicator_sum - len(second_state_set)

    return inclusion_indicator


def add_path_to_model(G, model, path_len, first_state_vars, last_state_vars, v_funcs=None):
    """
    Adds a path from first_state_vars to last_state_vars to the model, i.e. requires that last_state_vars
    represents the state resulting after path_len time steps from first_state_vars
    :param G:
    :param model:
    :param path_len:
    :param first_state_vars:
    :param last_state_vars:
    :param v_funcs:
    :return:
    """
    start = time.time()
    n = len(G.vertices)
    assert path_len >= 1

    previous_state_vars = first_state_vars
    for l in range(path_len):
        next_state_vars = last_state_vars if l == path_len - 1 else [
            model.addVar(vtype=gurobipy.GRB.BINARY, name="transient_path_state_var_{}_{}".format(l, i))
            for i in range(n)]
        model.update()

        for i in range(n):
            v_func = v_funcs[i] if v_funcs is not None else G.vertices[i].function
            predecessor_indices = [u.index for u in G.vertices[i].predecessors()]
            predecessor_vars = [previous_state_vars[index] for index in predecessor_indices]
            if len(predecessor_vars) == 0:
                model.addConstr(previous_state_vars[i] == next_state_vars[i],
                                name="stable_constraint_{}_node_{}".format(l, i))
                model.update()
            else:
                add_truth_table_consistency_constraints(model, v_func, next_state_vars[i], predecessor_vars,
                                                        name_prefix="transient_path_step_{}vertex_{}".format(l, i))

    # print "Time taken to add path constraints:{:.2f} seconds".format(time.time() - start)


# noinspection PyArgumentList
def attractors_ilp_with_keys(G, max_len=None, max_num=None,
                             model_type_restriction=FunctionTypeRestriction.NONE, simplify_general_boolean=False,
                             slice_size=15):
    total_start = time.time()
    part_start = time.time()
    T = 2**len(G.vertices) if max_len is None else max_len
    P = 2**len(G.vertices) if max_num is None else max_num
    n = len(G.vertices)

    model = gurobipy.Model()
    a_matrix = numpy.matrix([[model.addVar(vtype=gurobipy.GRB.BINARY, name="a_{}_{}".format(p, t))
                              for t in range(T+1)] for p in range(P)])  # TODO: fill for t=-1, and simplify code
    v_matrix = numpy.array([[[model.addVar(vtype=gurobipy.GRB.BINARY, name="v_{}_{}_{}".format(i, p, t))
                              for t in range(T+1)] for p in range(P)] for i in range(n)])
    model.update()

    # print "Time taken for basic prep:{:.2f} seconds".format(time.time() - part_start)
    # part_start = time.time()

    for p, t in itertools.product(range(P), range(T + 1)):
        # assert activity var meaning (ACTIVITY_SWITCH, MONOTONE, IF_NON_ACTIVE)
        # for i in range(n):
        #     model.addConstr(a_matrix[p, t] >= v_matrix[i, p, t], name="activity_switch_{}_{}_{}".format(i, p, t))

        model.addConstr(n * a_matrix[p, t] >= sum(v_matrix[i, p, t] for i in range(n)),
                        name="activity_{}_{}".format(p, t))
        if t < T:
            model.addConstr(a_matrix[p, t + 1] >= a_matrix[p, t], name="monotone_{}_{}".format(p, t))
        # model.update()

    for p in range(P):
        # should be redundant now that we return a_matrix[p, T-1] as the optimization function.
        model.addConstr(a_matrix[p, T] <= a_matrix[p, T - 1], name="if_non_active_{}".format(p))
        if p != P - 1:
            model.addConstr(a_matrix[p, T] <= a_matrix[p + 1, T], name="active_order_{}".format(p))

    # print "Time taken for activity constraints preparation:{:.2f} seconds".format(time.time() - part_start)
    # part_start = time.time()

    for i in range(n):
        # assert stable
        in_degree = len(G.vertices[i].predecessors())
        if in_degree == 0:
            for p, t in itertools.product(range(P), range(T+1)):
                if G.vertices[i].function is None:
                    desired_val = v_matrix[i, p, t]
                else:
                    assert G.vertices[i].function in [False, True]
                    desired_val = int(G.vertices[i].function)
                    last_activity = a_matrix[p, t-1] if t > 0 else 0
                    # set initial state, for first active state
                    model.addConstr(v_matrix[i, p, t] <= desired_val + 1 - a_matrix[p, t] + last_activity,
                                    name="stable_constant_<=_{}_{}_{}".format(i, p, t))
                    model.addConstr(v_matrix[i, p, t] >= desired_val + a_matrix[p, t] - last_activity - 1,
                                    name="stable_constant_>=_{}_{}_{}".format(i, p, t))
                if t < T:
                    # if a[p,t+1] and a[p,t] then v_matrix[i,p,t+1]=v_matrix[i,p,t]
                    model.addConstr(v_matrix[i, p, t+1] <= desired_val + 2 - a_matrix[p, t + 1] - a_matrix[p, t],
                                    name="stable_<=_{}_{}_{}".format(i, p, t))
                    model.addConstr(v_matrix[i, p, t+1] >= desired_val + a_matrix[p, t + 1] + a_matrix[p, t] - 2,
                                    name="stable_>=_{}_{}_{}".format(i, p, t))
    state_keys = [[unique_state_keys([v_matrix[i, p, t] for i in range(n)], slice_size=slice_size) for t in range(T+1)]
                  for p in range(P)]
    predecessors_vars = numpy.array([[[[v_matrix[vertex.index, p, t] for vertex in G.vertices[i].predecessors()]
                                     for t in range(T+1)] for p in range(P)] for i in range(n)])
    vertices_f_vars_list = [None] * n
    for i in range(n):
        # assert consistent
        in_degree = len(G.vertices[i].predecessors())
        if in_degree != 0:
            find_model = G.vertices[i].function is None
            if (find_model and model_type_restriction == FunctionTypeRestriction.NONE) \
                or (not find_model and not isinstance(G.vertices[i].function,
                                                      logic.SymmetricThresholdFunction)):
                if find_model:
                    find_model_f_vars = []
                    for var_comb_index, var_combination in enumerate(
                            itertools.product((False, True), repeat=len(G.vertices[i].predecessors()))):
                        find_model_f_vars.append(model.addVar(vtype=gurobipy.GRB.BINARY, name="f_{}_{}".format(i, var_comb_index)))
                        model.update()
                else:
                    find_model_f_vars = None
                vertices_f_vars_list[i] = find_model_f_vars
                for p, t in itertools.product(range(P), range(T)):
                    if simplify_general_boolean and not find_model:
                        add_simplified_consistency_constraints(model, G.vertices[i].function, v_matrix[i, p, t + 1],
                                                               predecessors_vars[i, p, t],
                                                               "consistency_{}_{}_{}".format(i, p, t),
                                                               activity_variable=a_matrix[p, t])
                    else:
                        add_truth_table_consistency_constraints(model, G.vertices[i].function, v_matrix[i, p, t + 1],
                                                               predecessors_vars[i, p, t],
                                                               "consistency_{}_{}_{}".format(i, p, t),
                                                               activity_variable=a_matrix[p, t],
                                                               find_model_f_vars=find_model_f_vars)
            else:  # symmetric threshold / logic gate
                if find_model:
                    signs = []
                    for input_index in range(in_degree):
                        sign_var = model.addVar(vtype=gurobipy.GRB.BINARY, name="f_{}_signs_{}".format(i, input_index))
                        signs.append(sign_var)
                    if model_type_restriction == FunctionTypeRestriction.SYMMETRIC_THRESHOLD:
                        threshold_expression = model.addVar(vtype=gurobipy.GRB.INTEGER, lb=0, ub=in_degree + 1,
                                                            name="f_{}_threshold".format(i))
                    else:
                        threshold_indicator = model.addVar(
                            vtype=gurobipy.GRB.BINARY, name="f_{}_gate_indicator".format(i))
                        threshold_expression = 1 + threshold_indicator * (in_degree - 1)
                    model.update()
                else:
                    signs = G.vertices[i].function.signs
                    threshold_expression = G.vertices[i].function.threshold
                for p, t in itertools.product(range(P), range(T)):
                    input_sum_expression = 0
                    for predecessors_index, (sign, var) in enumerate(zip(signs, predecessors_vars[i, p, t])):
                        if find_model:
                            # signs are variables, so can't multiply. need to have a variable signed_input s.t.:
                            # signed_input = not(xor(sign, var))
                            z = model.addVar(vtype=gurobipy.GRB.BINARY, name="signed_input_var_{}_{}_{}_{}".
                                             format(i, p, t, predecessors_index))
                            model.update()
                            model.addConstr(z >= -sign - var + 1, name="signed_input_constr_type00_{}_{}_{}_{}".
                                            format(i, p, t, predecessors_index))
                            model.addConstr(z <= sign - var + 1, name="signed_input_constr_type01_{}_{}_{}_{}".
                                            format(i, p, t, predecessors_index))
                            model.addConstr(z <= -sign + var + 1, name="signed_input_constr_type10_{}_{}_{}_{}".
                                            format(i, p, t, predecessors_index))
                            model.addConstr(z >= sign + var - 1, name="signed_input_constr_type11_{}_{}_{}_{}".
                                            format(i, p, t, predecessors_index))
                            # model.update()
                            input_sum_expression += z
                        else:
                            input_sum_expression += var if sign else (1 - var)
                    model.addConstr((in_degree + 1)*(v_matrix[i, p, t + 1] + 1 - a_matrix[p, t]) >=
                                    (input_sum_expression - threshold_expression + 1),
                                    name="monotone_func_consistency_>=_{}_{}_{}".format(i, p, t))
                    model.addConstr((in_degree + 1)*(v_matrix[i, p, t + 1] - 1 + a_matrix[p, t]) <=
                                    (in_degree + input_sum_expression - threshold_expression + 1),
                                    name="monotone_func_consistency_<=_{}_{}_{}".format(i, p, t))

    # print "Time taken for consistent and stable preparation:{:.2f} seconds".format(time.time() - part_start)
    # part_start = time.time()

    # CYCLIC
    for i, p, t in itertools.product(range(n), range(P), range(T)):
        # (~a[p, t - 1] & a[p, t]) >> EQ(p, p, t, T)
        last_activity = a_matrix[p, t - 1] if t > 0 else 0
        model.addConstr(a_matrix[p, t] - last_activity - 1 + v_matrix[i, p, t] <= v_matrix[i, p, T],
                        name="cyclic_<=_{}_{}_{}".format(i, p, t))
        model.addConstr(-a_matrix[p, t] + last_activity + 1 + v_matrix[i, p, t] >= v_matrix[i, p, T],
                        name="cyclic_>=_{}_{}_{}".format(i, p, t))

    # model.update()

    # print "Time taken for cyclic constraints preparation:{:.2f} seconds".format(time.time() - part_start)
    # part_start = time.time()

    # SIMPLE
    for t, p in itertools.product(range(T-1), range(P)):
        # (a[p, t] & a[p, t-1]) >> ~EQ(p, p, t, T) (assumes a_matrix monotone)
        strictly_larger_ind = create_state_keys_comparison_var(model=model, first_state_keys=state_keys[p][T-1],
                                                               second_state_keys=state_keys[p][t],
                                                               include_equality=False,
                                                               upper_bound=2**slice_size,
                                                               name_prefix="simple>_{}_{}".format(p, t))
        model.addConstr(strictly_larger_ind >= a_matrix[p, t],
                        name="simple_{}_{}".format(p, t))
    # model.update()

    # print "Time taken for simple constraints preparation:{:.2f} seconds".format(time.time() - part_start)
    # part_start = time.time()

    # UNIQUE
    # To reduce symmetry, using the order defined by a unique state id function,
    # constraint each attractor to have its final state be the largest one, and
    # constraint the final states of attractors to be monotone increasing amongst.
    # Also requires the k active attractors to be the last ones. TODO: verify
    # SHOULD result in only one optimal solution.
    for p in range(P-1):
        # as long as keys are uniques, this forces uniqueness if p and p + 1 are both active
        strictly_larger_ind = create_state_keys_comparison_var(model=model,
                                                               first_state_keys=state_keys[p + 1][T-1],
                                                               second_state_keys=state_keys[p][T-1],
                                                               include_equality=False,
                                                               upper_bound=2**slice_size,
                                                               name_prefix="key_order_between_attractors_{}".
                                                               format(p))
        # print "state keys for p={}".format(p)
        # print strictly_larger_ind.VarName
        # for key in state_keys[p+1][T]:
        #     print key
        # print "\n"
        # for key in state_keys[p][T]:
        #     print key
        # print "\n\n"
        model.addConstr(strictly_larger_ind >= a_matrix[p, T-1],  # need only a[p, T] because they're monotone
                        name="key_order_between_attractors_{}".format(p))
        # model.update()

    # Constraint the number of active attractors using 2**#input_nodes, P as lower and upper bounds.
    # lower bound can only be used if the maximal theoretical attractor length is allowed.
    model.addConstr(sum(a_matrix[p, T] for p in range(P)) <= P, name="upper_objective_bound")
    model.update()

    # print "Time taken for unique constraints preparation:{:.2f} seconds".format(time.time() - part_start)

    # print_model_constraints(model)
    # print model
    # print "Time taken for model preparation:{:.2f} seconds".format(time.time() - total_start)
    return model, a_matrix, v_matrix, state_keys, vertices_f_vars_list


def bitchange_attractor_ilp_with_keys(G, max_len=None, slice_size=15):
    T = max_len
    n = len(G.vertices)

    model = gurobipy.Model()
    a_list = [model.addVar(vtype=gurobipy.GRB.BINARY, name="a_{}".format(t)) for t in range(T+1)]
    v_matrix = numpy.array([[model.addVar(vtype=gurobipy.GRB.BINARY, name="v_{}_{}".format(i, t))
                              for t in range(T+1)] for i in range(n)])
    model.update()

    for t in range(T + 1):
        # assert activity var meaning (ACTIVITY_SWITCH, MONOTONE, IF_NON_ACTIVE)

        model.addConstr(n * a_list[t] >= sum(v_matrix[i, t] for i in range(n)),
                        name="activity_{}".format(t))
        if t < T:
            model.addConstr(a_list[t + 1] >= a_list[t], name="monotone_{}".format(t))

    model.addConstr(a_list[T] <= a_list[T - 1], name="if_non_active")
    model.addConstr(a_list[T] == 1, name="must_be_active")

    for i in range(n):
        # assert stable
        in_degree = len(G.vertices[i].predecessors())
        if in_degree == 0:
            for t in range(T+1):
                if G.vertices[i].function is None:
                    desired_val = v_matrix[i, t]
                else:
                    assert G.vertices[i].function in [False, True]
                    desired_val = int(G.vertices[i].function)
                    last_activity = a_list[t-1] if t > 0 else 0
                    # set initial state, for first active state
                    model.addConstr(v_matrix[i, t] <= desired_val + 1 - a_list[t] + last_activity,
                                    name="stable_constant_<=_{}_{}".format(i, t))
                    model.addConstr(v_matrix[i, t] >= desired_val + a_list[t] - last_activity - 1,
                                    name="stable_constant_>=_{}_{}".format(i, t))
                if t < T:
                    # if a[p,t+1] and a[p,t] then v_matrix[i,p,t+1]=v_matrix[i,p,t]
                    model.addConstr(v_matrix[i, t+1] <= desired_val + 2 - a_list[t + 1] - a_list[t],
                                    name="stable_<=_{}_{}".format(i, t))
                    model.addConstr(v_matrix[i, t+1] >= desired_val + a_list[t + 1] + a_list[t] - 2,
                                    name="stable_>=_{}_{}".format(i, t))
    state_keys = [unique_state_keys([v_matrix[i, t] for i in range(n)], slice_size=slice_size) for t in range(T+1)]
    predecessors_vars = numpy.array([[[v_matrix[vertex.index, t] for vertex in G.vertices[i].predecessors()]
                                     for t in range(T+1)] for i in range(n)])

    all_functions_bitchange_vars = []

    for i in range(n):
        # assert consistent
        in_degree = len(G.vertices[i].predecessors())
        if in_degree != 0:
            vertex_bitchange_vars = []
            find_model_f_vars = []
            for var_comb_index, var_combination in enumerate(
                    itertools.product((False, True), repeat=len(G.vertices[i].predecessors()))):
                bitchange_var = model.addVar(vtype=gurobipy.GRB.BINARY, name="f_{}_{}".format(i, var_comb_index))
                model.update()
                vertex_bitchange_vars.append(bitchange_var)
                all_functions_bitchange_vars.append(bitchange_var)
                find_model_f_vars.append(1 - bitchange_var if G.vertices[i].function(*var_combination) else bitchange_var)
            for t in range(T):
                add_truth_table_consistency_constraints(model, G.vertices[i].function, v_matrix[i, t + 1],
                                                       predecessors_vars[i, t],
                                                       "consistency_{}_{}".format(i, t),
                                                       activity_variable=a_list[t],
                                                       find_model_f_vars=find_model_f_vars)

    # CYCLIC
    for i, t in itertools.product(range(n), range(T)):
        # (~a[t - 1] & a[t]) >> EQ(t, T)
        last_activity = a_list[t - 1] if t > 0 else 0
        model.addConstr(a_list[t] - last_activity - 1 + v_matrix[i, t] <= v_matrix[i, T],
                        name="cyclic_<=_{}_{}".format(i, t))
        model.addConstr(-a_list[t] + last_activity + 1 + v_matrix[i, t] >= v_matrix[i, T],
                        name="cyclic_>=_{}_{}".format(i, t))

    # SIMPLE
    for t in range(T-1):
        # (a[t] & a[t-1]) >> ~EQ(t, T) (assumes a_list monotone)
        strictly_larger_ind = create_state_keys_comparison_var(model=model, first_state_keys=state_keys[T-1],
                                                               second_state_keys=state_keys[t],
                                                               include_equality=False,
                                                               upper_bound=2**slice_size,
                                                               name_prefix="simple>_{}".format(t))
        model.addConstr(strictly_larger_ind >= a_list[t],
                        name="simple_{}".format(t))

    model.update()

    # print_model_constraints(model)
    # print model
    return model, numpy.array([state_keys]), numpy.array([a_list]), all_functions_bitchange_vars


def set_mip_start(model, v_matrix, final_states_a_vars, attractors):
    """
    Given a model with its variable matrix and a set of pre-computed attractors, sets a MIP start
    to the model with those attractors as the first ones, and all others shut off.
    :param model:
    :param v_matrix: v_matrix[i, p, t] is activity of node i in attractor p at time t.
    :param attractors: A (non-consumable) iterable of attractors, each an iterable of states.
    Attractors are assumed to be bounded in length by the same bound in the model, and be less than the
    number of possible model attractors.
    :return:
    """
    n = v_matrix.shape[0]
    P = v_matrix.shape[1]
    T = v_matrix.shape[2] - 1
    order_key_func = lambda node_states: sum(node * 2**i for (i, node) in enumerate(node_states))

    # Rotate the attractors to meet the order used in the ILP - last state should be largest.
    ordered_attractors = list()
    assert len(attractors) <= P
    for p in range(len(attractors)):
        length = len(attractors[p])
        assert length <= T
        largest_ind = max(range(length), key=lambda t: order_key_func(attractors[p][t]))
        ordered_attractors.append([attractors[p][(t + largest_ind + 1) % length] for t in range(length)])

    # Sort the attractor list in increasing order.
    ordered_attractors.sort(key=lambda attractor: order_key_func(attractor[-1]), reverse=False)

    # set model values to undefined before setting the defined ones (redundant?)
    for var in model.getVars():
        var.start = GRB.UNDEFINED

    # Assign attractors to last columns of v_matrix.
    for p in range(P):
        sim_p = p - (P - len(ordered_attractors))

        if sim_p < 0:
            for i, t in itertools.product(range(n), range(T + 1)):
                v_matrix[i, p, t].start = 0
            final_states_a_vars[p].start = 0
        else:
            length = len(ordered_attractors[sim_p])
            for i, t in itertools.product(range(n), range(T - length)):
                v_matrix[i, p, t].start = 0
            for i, t in itertools.product(range(n), range(T - length, T)):
                sim_t = t - (T - length)
                v_matrix[i, p, t].start = ordered_attractors[sim_p][sim_t][i]
            for i in range(n):
                v_matrix[i, p, T].start = ordered_attractors[sim_p][0][i]
            final_states_a_vars[p].start = 1
    model.Params.StartNodeLimit = 2000000000 # maximal possible value
    # model.Params.Heuristics = 0 # TODO: remove after I have the ignoring mip start issue solved.
    model.update()
    # v_string = ""
    # for p in range(P):
    #     v_string += "p={}\n".format(p)
    #     for t in range(T + 1):
    #         for i in range(n):
    #             v_string += str(v_matrix[i, p, t].start) + ", "
    #         v_string += "\t t={}\n".format(t)
    # print v_string
    return None


def get_expr_coos(expr, var_indices):
    for i in range(expr.size()):
        dvar = expr.getVar(i)
        yield expr.getCoeff(i), var_indices[dvar]


def get_matrix_coos(m):
    dvars = m.getVars()
    constrs = m.getConstrs()
    var_indices = {v: i for i, v in enumerate(dvars)}
    indices_to_vars = {i: v for (v, i) in var_indices.items()}
    for row_idx, constr in enumerate(constrs):
        for coeff, col_idx in get_expr_coos(m.getRow(constr), var_indices):
            yield row_idx, constr.ConstrName, indices_to_vars[col_idx].VarName, coeff


def print_model_values(model, model_vars=None):
    # assumes optimization has completed successfully
    if not model_vars:
        model_vars = model.getVars()
    for var in model_vars:
        print "{}\t{}".format(var.VarName, var.X)


def print_attractors(model):
    """
    Prints the attractors of a model after it has been optimized.
    :param model:
    :return:
    """
    v_name_parts = [var.VarName.split("_")[1:] for var in model.getVars() if "v_" in var.VarName]
    max_i = 0
    max_p = 0
    max_t = 0
    for name_part_list in v_name_parts:
        assert len(name_part_list) == 3
        i, p, t = name_part_list
        max_i = max(max_i, int(i))
        max_p = max(max_p, int(p))
        max_t = max(max_t, int(t))
    T = max_t
    P = max_p + 1
    n = max_i + 1
    a_variables = [[model.getVarByName("a_{}_{}".format(p, t)) for t in range(T+1)] for p in range(P)]
    v_variables = [[[model.getVarByName("v_{}_{}_{}".format(i, p, t)) for t in range(T)] for p in range(P)]
                   for i in range(n)]
    # assert len([None for p in range(P) if int(round(a_variables[p][T-1].X)) == 1]) == n_attractors
    n_attractors = sum(last_var.X for last_var in [a_list[-1] for a_list in a_variables])
    for attractor_number in range(int(n_attractors)):
        p = P - attractor_number - 1
        length = len([None for t in range(T) if int(round(a_variables[p][t].X)) == 1])
        print "Attractor #{}, length {}".format(attractor_number + 1, length)
        # TODO: support for original graph names?
        # print reduce(lambda a, b: "{}\t{}".format(a, b), ["v_{}".format(i) for i in range(n)])
        for t in range(T - length, T):
            # noinspection PyTypeChecker
            print reduce(lambda a, b: "{}{}".format(a, b), [int(round(v_variables[i][p][t].X)) for i in range(n)])


def print_attractors_enumeration(model):
    """
    Prints the attractors of a model after it has been optimized, when
    each solution corresponds to one attractor.
    :param model:
    :return:
    """
    n_attracotrs = model.SolCount
    v_name_parts = [var.VarName.split("_")[1:] for var in model.getVars() if "v_" in var.VarName]
    max_i = 0
    max_t = 0
    for name_part_list in v_name_parts:
        assert len(name_part_list) == 3
        i, p, t = name_part_list
        max_i = max(max_i, int(i))
        max_t = max(max_t, int(t))
    T = max_t
    P = 1
    n = max_i + 1
    a_variables = [[model.getVarByName("a_{}_{}".format(p, t)) for t in range(T+1)] for p in range(P)]
    v_variables = [[[model.getVarByName("v_{}_{}_{}".format(i, p, t)) for t in range(T)] for p in range(P)]
                   for i in range(n)]
    for p in range(n_attracotrs):
        model.setParam(gurobipy.GRB.Param.SolutionNumber, p)
        length = len([None for t in range(T) if int(round(a_variables[0][t].Xn)) == 1])
        print "Attractor #{}, length {}".format(p + 1, length)
        # TODO: support for original graph names?
        # print reduce(lambda a, b: "{}\t{}".format(a, b), ["v_{}".format(i) for i in range(n)])
        for t in range(T - length, T):
            # noinspection PyTypeChecker
            print reduce(lambda a, b: "{}{}".format(a, b), [int(round(v_variables[i][0][t].Xn)) for i in range(n)])


def get_model_attractors(model):
    """
    Returns the attractors of a model, after it has been optimized, when
    each solution corresponds to one attractor.
    :param model:
    :return:
    """
    n_attracotrs = model.SolCount
    v_name_parts = [var.VarName.split("_")[1:] for var in model.getVars() if "v_" in var.VarName]
    max_i = 0
    max_t = 0
    for name_part_list in v_name_parts:
        assert len(name_part_list) == 3
        i, p, t = name_part_list
        max_i = max(max_i, int(i))
        max_t = max(max_t, int(t))
    T = max_t
    P = 1
    n = max_i + 1
    a_variables = [[model.getVarByName("a_{}_{}".format(p, t)) for t in range(T+1)] for p in range(P)]
    v_variables = [[[model.getVarByName("v_{}_{}_{}".format(i, p, t)) for t in range(T)] for p in range(P)]
                   for i in range(n)]
    attractors = []
    for p in range(n_attracotrs):
        attractor = []
        model.setParam(gurobipy.GRB.Param.SolutionNumber, p)
        length = len([None for t in range(T) if int(round(a_variables[0][t].Xn)) == 1])
        for t in range(T - length, T):
            state = [int(round(v_variables[i][0][t].Xn)) for i in range(n)]
            attractor.append(state)
        attractors.append(attractor)
    return attractors


def print_model_constraints(model):
    quadruples = list(get_matrix_coos(model))
    constr_attrs = [(constr.Sense, constr.RHS) for constr in model.getConstrs()]
    for constraint_index in range(max(row_index for row_index, _, _, _ in quadruples) + 1):
        con_str = None
        for row_index, constr_name, var_name, coeff in quadruples:
            if row_index == constraint_index:
                if not con_str:
                    con_str = constr_name + ": "
                con_str += " {}{}{}".format(str(coeff) if coeff not in [1.0, -1.0] else "",
                                            "-" if coeff == -1.0 else "+" if coeff == 1 else "", var_name)
        con_str += " {} {}".format(*constr_attrs[constraint_index])
        print con_str


def print_opt_solution(model):
    name_val_pairs = []
    for var in model.getVars():
        name_val_pairs.append((var.VarName, var.X))
    val_str = ""
    for pair in name_val_pairs:
        val_str += "{} = {}\n".format(pair[0], int(pair[1]))
    print val_str


def steady_state_ilp(G, simplify_general_boolean=False):
    model = gurobipy.Model()
    v_dict = dict()
    for v in G.vertices:
        if len(v.predecessors()) == 0 and v.function is not None:
            v_dict[v] = int(bool(v.function))
        else:
            v_dict[v] = model.addVar(vtype=gurobipy.GRB.BINARY, name=v.name)
    model.update()

    for i, v in enumerate(G.vertices):
        if len(v.predecessors()) == 0 or isinstance(v_dict[v], int):
            continue
        func = v.function
        predecessor_vars = [v_dict[u] for u in v.predecessors()]
        if simplify_general_boolean:
            add_simplified_consistency_constraints(model, v.function, v_dict[v], predecessor_vars,
                                                   "consistency_{}".format(v.name))
        else:  # TODO: support find model and other stuff
            add_truth_table_consistency_constraints(model, v.function, v_dict[v], predecessor_vars,
                                                    "consistency_{}".format(v.name))
    return model, v_dict
