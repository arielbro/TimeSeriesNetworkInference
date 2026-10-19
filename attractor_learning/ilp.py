import gurobipy
import itertools
import time
import math
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


def _binary_constant(value):
    value = int(round(float(value)))
    if value not in (0, 1):
        raise ValueError("v must be binary, but got value {}".format(value))
    return value


def _agreeing_input(model, positive, negative, value, name):
    """1 iff an input's value agrees with its edge's sign: the activating-edge binary p for an input at 1, the
    inhibiting-edge binary q for one at 0. For a variable input v this is p*v + q*(1 - v), written as its convex
    hull given p + q <= 1 (one binary, five rows)."""
    if not isinstance(value, gurobipy.Var):
        return positive if _binary_constant(value) == 1 else negative
    agrees = model.addVar(vtype=GRB.BINARY, name=name)
    model.addConstr(agrees <= value + negative, name=name + "_hull1")
    model.addConstr(agrees <= 1 - value + positive, name=name + "_hull2")
    model.addConstr(agrees <= positive + negative, name=name + "_hull3")
    model.addConstr(agrees >= positive + value - 1, name=name + "_hull4")
    model.addConstr(agrees >= negative - value, name=name + "_hull5")
    return agrees


def add_state_agreement_terms(model, predicted_state, target_state, prefix=""):
    """Per-cell agreement of a predicted state (variables) with a target state, as terms for an objective that
    maximizes their sum. Against a constant target cell the term is the predicted cell itself (or its
    complement). Against a variable one (a denoised data row) it is a binary z bounded by z <= 1 + x - y and
    z <= 1 + y - x: the objective maximizes it, so the two rows bounding it from below are not needed - which
    is also why these terms are only valid in a maximized objective, with a positive weight."""
    terms = []
    for index, (predicted, target) in enumerate(zip(predicted_state, target_state)):
        if not isinstance(target, gurobipy.Var):
            terms.append(predicted if _binary_constant(target) == 1 else 1 - predicted)
            continue
        name = "{}_agreement_index_{}".format(prefix, index)
        agrees = model.addVar(vtype=GRB.BINARY, name=name)
        model.addConstr(agrees <= 1 + predicted - target, name=name + "_<=1")
        model.addConstr(agrees <= 1 + target - predicted, name=name + "_<=2")
        terms.append(agrees)
    return terms


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


def add_path_to_model(G, model, path_len, first_state_vars, model_f_vars, max_indegree=None, name_prefix=""):
    """
    Adds a path of path_len synchronous steps from first_state_vars (variables or constant 0/1 values) under
    signed symmetric threshold functions, and returns the states along it as new binaries: excluding the first
    state, including the last.
    model_f_vars gives each node's function as (positives, negatives, threshold, unused), aligned with the node's
    predecessors in G: per input, a binary that is 1 iff the input is used and activating and one that is 1 iff
    it is used and inhibiting (at most one of which may be 1, which the caller constrains), the threshold, and
    an indicator that is 1 iff the node uses no input at all, in which case it holds its value, as an input node
    does in Network.next_state. unused may be None for a node that always uses its inputs. Any of them may be
    constants, to check a fixed model. A node without predecessors holds its value.
    :param max_indegree: upper bound on how many inputs a node may end up using, when the caller constrains it.
        Only tightens the big-M of the threshold comparison, which is otherwise sized by the predecessor count.
    """
    assert path_len >= 1, "can't add path constraint with path of length {}".format(path_len)
    n = len(G.vertices)
    previous_state = first_state_vars
    states = []
    for step in range(path_len):
        prefix = "{}_step_{}".format(name_prefix, step)
        next_state = [model.addVar(vtype=GRB.BINARY, name="{}_state_{}".format(prefix, i)) for i in range(n)]
        for i in range(n):
            predecessor_indices = [u.index for u in G.vertices[i].predecessors()]
            if not predecessor_indices:
                model.addConstr(next_state[i] == previous_state[i], name="{}_stable_{}".format(prefix, i))
                continue
            positives, negatives, threshold, unused = model_f_vars[i]
            agreeing = [_agreeing_input(model, p, q, previous_state[k], "{}_agree_{}_{}".format(prefix, i, k))
                        for k, p, q in zip(predecessor_indices, positives, negatives)]
            # strictly positive iff the output is 1. It lies in [2 - degree, degree] for a used node and is
            # exactly 1 for an unused one (no inputs, threshold 0), so the >= side only has to be relaxed by
            # unused and the <= side not at all; the indicator below then holds the node's value instead.
            comparison = gurobipy.quicksum(agreeing) - threshold + 1
            # the degree cannot exceed the predecessor count nor, where the caller imposes one, max_indegree
            degree_bound = len(predecessor_indices) if max_indegree is None \
                else min(len(predecessor_indices), max_indegree)
            big_m = degree_bound + 1
            slack = 0 if unused is None else unused
            model.addConstr(big_m * next_state[i] >= comparison - slack,
                            name="{}_threshold_function_path_constraint_>=_{}".format(prefix, i))
            model.addConstr(big_m * next_state[i] <= big_m + comparison - 1,
                            name="{}_threshold_function_path_constraint_<=_{}".format(prefix, i))
            if unused is not None:
                previous = previous_state[i] if isinstance(previous_state[i], gurobipy.Var) \
                    else _binary_constant(previous_state[i])
                model.addGenConstrIndicator(unused, True, next_state[i] == previous,
                                            name="{}_unused_node_holds_value_{}".format(prefix, i))
        states.append(next_state)
        previous_state = next_state
    return states


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
