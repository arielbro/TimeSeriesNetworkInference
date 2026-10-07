import gurobipy
import itertools
import time
import math
from attractor_learning.graphs import FunctionTypeRestriction
from gurobipy import GRB

# TODO: find good upper bound again, why didn't 29 work on MAPK_large2?
# http://files.gurobi.com/Numerics.pdf a good resource on numerical issues, high values cause them.
# TODO: gurobipy now supports adding or, and, if-then constraints (among others). I should rewrite what I implemented by hand to use those interfaces.


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

    # print("parts={}".format(parts))
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
    assert len(first_state_keys) == len(second_state_keys), "length of state key lists differ"
    last_var = 0 if not include_equality else 1
    M = upper_bound  # M - 1 is actual largest value
    for i in range(len(first_state_keys)):
        a = first_state_keys[-i - 1]
        b = second_state_keys[-i - 1]
        z = model.addVar(vtype=gurobipy.GRB.BINARY, name="{}_{}_indicator".format(name_prefix, i))
        model.update()
        model.addConstr(M * z >= a - b + last_var, name="{}_{}_>=constraint".format(name_prefix, i))
        model.addConstr(M * z <= a - b + last_var + (M - 1), name="{}_{}_<=constraint".format(name_prefix, i))
        last_var = z
        # print("a_{}={}, b_{}={}, M={}, M+1={}".format(len(first_state_keys) - i-1, a, len(first_state_keys) -i-1, b, M, M+1))
    return last_var


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
        elif isinstance(v_func(*var_combination), gurobipy.Var):
            desired_val = v_func(*var_combination)
        else:
            desired_val = 1 if v_func(*var_combination) else 0  # == because sympy
        # this expression is |in_degree| iff their states agrees with var_combination
        indicator_expression = gurobipy.quicksum(v if state else 1 - v for (v, state) in
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


def add_state_equality_indicator(model, first_state, second_state,
                                 force_equal=True, prefix=None, return_per_index=False):
    """
    Adds and returns a (NON STRICT) binary indicator for whether two network states are equal.
    result = 1 -> first_state == second_state.
    If force_equal, result = 1 <-> first_state == second_state.
    Does not use state hashing and ordering. States are allowed to be variables or constants.
    If return_per_index, returns the list of per-node (per-cell) equality indicators instead of a
    single whole-state indicator, allowing cell-wise rather than row-wise agreement scoring.

    :param model:
    :param first_state:
    :param second_state:
    :param force_equal:
    :param prefix:
    :param return_per_index:
    :return:
    """
    # TODO: write tests
    equality_indicators = []
    for index, (first_val, second_val) in enumerate(zip(first_state, second_state)):
        # add indicator for whether these two values are equal.
        val_equality_indicator = model.addVar(vtype=GRB.BINARY, name="{}_state_equality_indicator_index_{}".
                                                format(prefix, index))
        if force_equal:
            # force eq_var = 1 <-> first_val == second_val
            model.addConstr(val_equality_indicator <= 1 + second_val - first_val,
                                name="{}_state_equality_constraint_<=_1_index_{}".
                                format(prefix, index))
            model.addConstr(val_equality_indicator <= 1 + first_val - second_val,
                                name="{}_state_equality_constraint_<=_2_index_{}".
                                format(prefix, index))
            model.addConstr(val_equality_indicator >= first_val + second_val - 1,
                                name="{}_state_equality_constraint_>=_1_index_{}".
                                format(prefix, index))
            model.addConstr(val_equality_indicator >= 1 - first_val - second_val,
                                name="{}_state_equality_constraint_>=_2_index_{}".
                                format(prefix, index))
        else:
            model.addGenConstrIndicator(val_equality_indicator, 1, first_val == second_val,
                                        name="{}_state_equality_constraint_index_{}".format(prefix, index))
        equality_indicators.append(val_equality_indicator)
    if return_per_index:
        return equality_indicators
    equality_var = model.addVar(vtype=GRB.BINARY, name="{}_state_equality_indicator".
                 format(prefix))
    model.update()
    model.addGenConstrAnd(equality_var, equality_indicators,
                          name="{}_state_inequality_constraint".format(prefix))
    model.update()
    return equality_var


def add_hashed_state_inclusion_indicator(model, first_state, second_state_set, slice_size, prefix=None,
                                         assume_uniqueness=True):
    """
    Adds and returns a binary indicator for whether one network state is included in a set of others.
    States should be either iterables of constant binary values (False/True/0/1), or model variables.
    If assume_uniqueness, it is assumed that all states in second_state_set are unique, and might not
    return a binary indicator otherwise.
    Indicator is created with use of state hashing and order indicators for state pairs.
    :param model:
    :param first_state:
    :param second_state_set:
    :return:
    """
    indicator_sum = 0
    first_state_keys = unique_state_keys(first_state, slice_size=slice_size)
    # print("len of second state set - {}".format(len(second_state_set)))
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
        # print("indicator sum pre: {}".format(indicator_sum))
        indicator_sum += larger_var + smaller_var
        # print("indicator sum post: {}".format(indicator_sum))
    model.update()

    # We want an indicator for inclusion of first_state in the second state, which is equivalent to its equality
    # with at least one state there (exactly one if they are unique), or at least len(second_state_set) + 1
    # indicators with value 1.
    diff = indicator_sum - len(second_state_set)
    if assume_uniqueness:
        # since each state admits one indicator with value 1, and at most one admits 2. So the difference is binary.
        inclusion_indicator = diff
    else:
        # we want an indicator for whether diff > 0. Note that diff is in range [0, len(second_state_set)]
        inclusion_indicator = model.addVar(vtype=gurobipy.GRB.BINARY,
                                           name="{}_inclusion_indicator".format(prefix))
        model.update()
        model.addConstr(diff >= inclusion_indicator,
                        name="{}_inclusion_indicator_constraint_>=".format(prefix))
        model.addConstr(diff <= inclusion_indicator * len(second_state_set),
                        name="{}_inclusion_indicator_constraint_<=".format(prefix))

    return inclusion_indicator


def add_path_to_model(G, model, path_len, first_state_vars, model_f_vars, last_state_vars=None,
                      v_funcs_restrictions=None, name_prefix="", max_indegree=None):
    """
    Adds a path from first_state_vars to last_state_vars to the model, i.e. requires that last_state_vars
    represents the state resulting after path_len time steps from first_state_vars.
    Returns the state variables representing the path, excluding the first state and including the last state.
    If last_state_vars is undefined, creates a new state for it, and return it with the rest.
    If v_funcs_restrictions is not None, assumes it gives each function a possible restriction to symmetric threshold
    or a simple gate, with a given representation in model_f_vars.
    :param G:
    :param model:
    :param path_len:
    :param first_state_vars:
    :param last_state_vars:
    :param model_f_vars:
    :param max_indegree: upper bound on how many inputs a node may end up using, when the caller constrains
        it (the symmetric-threshold formulation lets a node zero out inputs it doesn't use). Only tightens
        the big-M of the threshold comparison, which is otherwise sized by the candidate count; None leaves
        it at the candidate count.
    :return:
    """
    start = time.time()
    n = len(G.vertices)
    assert path_len >= 1, "can't add path constraint with path of length {}".format(path_len)

    previous_state_vars = first_state_vars
    new_state_vars_list = []
    for l in range(path_len):
        next_state_vars = last_state_vars if ((last_state_vars is not None) and (l == path_len - 1)) else [
            model.addVar(vtype=gurobipy.GRB.BINARY, name="transient_path_state_var_{}_{}".format(l, i))
            for i in range(n)]
        new_state_vars_list.append(next_state_vars)

        for i in range(n):
            if len(G.vertices[i].predecessors()) == 0:
                model.addConstr(previous_state_vars[i] == next_state_vars[i],
                                name="stable_constraint_{}_node_{}".format(l, i))
            else:
                predecessor_indices = [u.index for u in G.vertices[i].predecessors()]
                predecessor_vars = [previous_state_vars[index] for index in predecessor_indices]

                if (v_funcs_restrictions is not None) and (
                        v_funcs_restrictions[i] == FunctionTypeRestriction.SYMMETRIC_THRESHOLD):
                    signs, threshold = model_f_vars[i][0], model_f_vars[i][1]
                    # optional third element: an indicator that is 1 iff the node ends up using no input at
                    # all (all its signs are zero). Such a node holds its value, exactly as an input node
                    # does in Network.next_state, instead of being constant True (an empty signed sum meets
                    # the threshold, which is forced to 0 at degree 0).
                    unused_node = model_f_vars[i][2] if len(model_f_vars[i]) > 2 else None
                    # v * s + (1-s)/2 gives v if s=1 and (1-v) if s=-1
                    # !!! This only works if s\in {-1, 1}, if we're learning the topology and s=0 is possible, breaks.
                    # instead, divide to cases based on the value of v and an extra variable.
                    # TODO: consolidate formulation here and in other parts of ilp.py
                    signed_input_vars = [model.addVar(vtype=gurobipy.GRB.BINARY,
                                         name="{}_signed_input_var_{}_{}_{}".format(name_prefix, l, i, j))
                                         for j in range(len(predecessor_vars))]
                    for j, (sign, var, signed_input) in enumerate(zip(signs, predecessor_vars, signed_input_vars)):
                        # we want signed_input = 1 <-> (s = 1 and v = 1) or (s = -1 and v = 0)
                        # v is binary, and s can take -1, 1
                        # (2v - 1)s yields 1 when v and s "agree" (1 and 1 or 0 and -1), otherwise less.
                        # And so, forcing signed_input to be 1 iff (2v - 1)s is 1 works, but that'd probably be harder
                        # to solve (since constraint-heavy expressions generally yield better relaxations).
                        # instead, just split to cases, since v is a constant
                        if isinstance(var, gurobipy.Var):
                            # v is a decision variable (e.g. a flippable/denoised data cell), so we can't
                            # split on its value. Linearize the product prod = sign * v via McCormick
                            # envelopes (exact since v is binary and sign in [-1, 1]):
                            #   v=1 -> prod=sign ; v=0 -> prod=0.
                            # Then (2v - 1) * sign == 2 * prod - sign, which lies in {-1, 0, 1} and equals
                            # 1 exactly when v and sign "agree", so signed_input is its positive part.
                            prod = model.addVar(lb=-1, ub=1, vtype=gurobipy.GRB.INTEGER,
                                                name="{}_signed_input_prod_{}_{}_{}".format(name_prefix, l, i, j))
                            model.addConstr(prod <= var, name="{}_signed_input_prod_mc1_{}_{}_{}".format(name_prefix, l, i, j))
                            model.addConstr(prod >= -var, name="{}_signed_input_prod_mc2_{}_{}_{}".format(name_prefix, l, i, j))
                            model.addConstr(prod <= sign + (1 - var), name="{}_signed_input_prod_mc3_{}_{}_{}".format(name_prefix, l, i, j))
                            model.addConstr(prod >= sign - (1 - var), name="{}_signed_input_prod_mc4_{}_{}_{}".format(name_prefix, l, i, j))
                            agreement = 2 * prod - sign  # in {-1, 0, 1}, == 1 iff v agrees with sign
                            model.addConstr(signed_input >= agreement,
                                            name="{}_signed_input_constraint_{}_{}_{}>=".format(name_prefix, l, i, j))
                            model.addConstr(2 * signed_input <= agreement + 1,
                                            name="{}_signed_input_constraint_{}_{}_{}<=".format(name_prefix, l, i, j))
                        elif var == 1:
                            # s = 1 -> signed_input = 1, otherwise signed_input = 0
                            model.addConstr(signed_input >= sign,
                                            name="{}_signed_input_constraint_{}_{}_{}>=".format(name_prefix, l, i, j))
                            model.addConstr(2 * signed_input <= 1 + sign,
                                            name="{}_signed_input_constraint_{}_{}_{}<=".format(name_prefix, l, i, j))
                        elif var == 0:
                            # s = -1 -> signed_input = 1, otherwise signed_input = 0
                            model.addConstr(signed_input >= -sign,
                                            name="{}_signed_input_constraint_{}_{}_{}>=".format(name_prefix, l, i, j))
                            model.addConstr(2 * signed_input <= 1 - sign,
                                            name="{}_signed_input_constraint_{}_{}_{}<=".format(name_prefix, l, i, j))
                        else:
                            raise ValueError("v must be binary, but got value {}".format(var))
                    # This expression ranges in [-threshold + 1, degree + 1]
                    # which is invariably bounded by [-len(predecessor_vars), len(predecessor_vars) + 1]
                    # It's strictly positive iff output should be 1.
                    input_to_threshold_comparison = gurobipy.quicksum(signed_input_vars) - threshold + 1
                    # big_m * unused_node relaxes both comparisons for a node that uses no input (the
                    # comparison then ranges in [-degree, degree + 1], so big_m makes them vacuous), leaving
                    # the indicator constraint below to hold its value instead.
                    # the comparison ranges within +/- (degree + 1), and the degree cannot exceed the
                    # candidate count nor, where the caller imposes one, max_indegree - so the smaller of
                    # the two is a valid and tighter big-M
                    degree_bound = len(predecessor_vars) if max_indegree is None                         else min(len(predecessor_vars), max_indegree)
                    big_m = degree_bound + 1
                    unused_slack = 0 if unused_node is None else big_m * unused_node
                    model.addConstr(big_m * next_state_vars[i] >= input_to_threshold_comparison - unused_slack,
                        name="{}_threshold_function_path_constraint_>=".format(name_prefix))
                    model.addConstr(big_m * next_state_vars[i] <=
                                    big_m + input_to_threshold_comparison - 1 + unused_slack,
                                    name="{}_threshold_function_path_constraint_<=".format(name_prefix))
                    if unused_node is not None:
                        model.addGenConstrIndicator(
                            unused_node, True, next_state_vars[i] == previous_state_vars[i],
                            name="{}_unused_node_holds_value_{}_{}".format(name_prefix, l, i))

                elif (v_funcs_restrictions is not None) and (
                        v_funcs_restrictions[i] == FunctionTypeRestriction.SIMPLE_GATES):
                    raise NotImplementedError()
                else:
                    assert (v_funcs_restrictions is None) or (v_funcs_restrictions[i] == FunctionTypeRestriction.NONE) \
                        or (v_funcs_restrictions[i] is None)
                    find_model_f_vars = model_f_vars[i] if (model_f_vars is not None) else None
                    v_func = None if ((model_f_vars is not None) and (model_f_vars[i]) is not None) \
                        else G.vertices[i].function
                    predecessor_indices = [u.index for u in G.vertices[i].predecessors()]
                    predecessor_vars = [previous_state_vars[index] for index in predecessor_indices]
                    add_truth_table_consistency_constraints(model, v_func, next_state_vars[i], predecessor_vars,
                                                            name_prefix="transient_path_step_{}vertex_{}".format(l, i),
                                                            activity_variable=None, find_model_f_vars=find_model_f_vars)

        previous_state_vars = next_state_vars

    # print("Time taken to add path constraints:{:.2f} seconds".format(time.time() - start))
    return new_state_vars_list


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
        print("{}\t{}".format(var.VarName, var.X))


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
        print(con_str)


def print_opt_solution(model):
    name_val_pairs = []
    for var in model.getVars():
        name_val_pairs.append((var.VarName, var.X))
    val_str = ""
    for pair in name_val_pairs:
        val_str += "{} = {}\n".format(pair[0], int(pair[1]))
    print(val_str)
