import gurobipy
import sympy
import numpy
import itertools
import logic
import time
import math
from graphs import FunctionTypeRestriction

unique_key_slicing_size = 10  # TODO: find elegant reformatting for this
# TODO: find good upper bound again, why didn't 29 work on MAPK_large2?
# http://files.gurobi.com/Numerics.pdf a good resource on numerical issues, high values cause them.
simplify_general_boolean = True # TODO: reformat as parameter, test both options (and benchmark times)


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


def unique_state_keys(ordered_state_variables):
    """
    Assign a unique numbering to a state, by summing 2**i over active vertices.
    Split the value among several variables, if needed, to fit in 32bit ints with room for mathematical operations.
    :param ordered_state_variables: an ordered fixed iterable of vertex state variables.
    :return: A key, a tuple of integers, identifying this state uniquely among all other states.
    """
    # TODO: see if possible to use larger slices.
    slice_size = unique_key_slicing_size
    # according to gurobi documentation (https://www.gurobi.com/documentation/7.5/refman/variables.html#subsubsection:IntVars),
    # int values are restrained to +-20 billion, or 2**30.9. I choose 2**29 as a very conservative bound.
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
    NOTE: creates len(first_state_keys) - 1 auxiliary variables. Assumes numbers are bounded by 2**30
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
    M = upper_bound
    for i in range(len(first_state_keys)):
        a = first_state_keys[-i - 1]
        b = second_state_keys[-i - 1]
        z = model.addVar(vtype=gurobipy.GRB.BINARY, name="{}_{}_indicator".format(name_prefix, i))
        model.update()
        model.addConstr((M + 1)*z >= a - b + last_var, name="{}_{}_>=constraint".format(name_prefix, i))
        model.addConstr((M + 1)*z <= a - b + last_var + M, name="{}_{}_<=constraint".format(name_prefix, i))
        last_var = z
        # print "a_{}={}, b_{}={}, M={}, M+1={}".format(len(first_state_keys) - i-1, a, len(first_state_keys) -i-1, b, M, M+1)
    return last_var


def add_truth_table_consistency_constraints(G, model, i, p, t, predecessors_vars,
                                            a_matrix, v_matrix, find_model_f_vars):
    """
    Adds consistency constraints to a model, as in Roded's paper.
    :param G:
    :param model:
    :param i:
    :param P:
    :param T:
    :param predecessors_vars:
    :param a_matrix:
    :param v_matrix:
    :param find_model:
    :return:
    """
    in_degree = len(G.vertices[i].predecessors())
    for var_comb_index, var_combination in enumerate(itertools.product((False, True), repeat=in_degree)):
        if find_model_f_vars:
            desired_val = find_model_f_vars[var_comb_index]
        else:
            desired_val = 1 if G.vertices[i].function(*var_combination) else 0  # == because sympy
        # this expression is |in_degree| iff their states agrees with var_combination
        indicator_expression = sum(v if state else 1 - v for (v, state) in
                                   zip(predecessors_vars[i, p, t], var_combination))
        # a[p, t] & (indicator_expression = in_degree) => v[i,p,t+1] = f(var_combination).
        # For x&y => a=b, require a <= b + (2 -x -b), a >= b - (2 -x -y)
        model.addConstr(v_matrix[i, p, t + 1] >= desired_val -
                        (in_degree + 1 - indicator_expression - a_matrix[p, t]),
                        name="consistent_>=_{}_{}_{}".format(i, p, t))
        model.addConstr(v_matrix[i, p, t + 1] <= desired_val +
                        (in_degree + 1 - indicator_expression - a_matrix[p, t]),
                        name="consistent_<=_{}_{}_{}".format(i, p, t))


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
        argument_vars = [build_logic_function_vars(argument, model, name_prefix + "_And_args_{}",
                                                   symbols_to_variables_dict)
                                                   for (i, argument) in enumerate(formula.args)]
        andVar = model.addVar(vtype=gurobipy.GRB.BINARY,
                              name="{}_And".format(name_prefix))
        for (i, argument_var) in enumerate(argument_vars):
            model.addConstr(andVar <= argument_var, name="{}_And_res_<=_{}".format(name_prefix, i))
        model.addConstr(andVar >= sum(argument_vars) - (len(argument_vars) - 1),
                        name="{}_And_res_>=".format(name_prefix, i))
        return andVar
    elif isinstance(formula, sympy.Or):
        argument_vars = [build_logic_function_vars(argument, model, name_prefix + "_Or_args_{}",
                                                   symbols_to_variables_dict)
                                                   for (i, argument) in enumerate(formula.args)]
        orVar = model.addVar(vtype=gurobipy.GRB.BINARY,
                              name="{}_Or".format(name_prefix))
        for (i, argument_var) in enumerate(argument_vars):
            model.addConstr(orVar >= argument_var, name="{}_Or_res_>=_{}".format(name_prefix, i))
        model.addConstr(orVar <= sum(argument_vars),
                        name="{}_And_res_<=".format(name_prefix, i))
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


def add_simplified_consistency_constraints(G, model, i, p, t, predecessors_vars,
                                               a_matrix, v_matrix):
    """
    Adds consistency constraints to the model by transforming each vertex' logic function to an indicator variable,
    and enforcing its value to be equal to the next state's value.
    :param G:
    :param model:
    :param i:
    :param P:
    :param T:
    :param predecessors_vars:
    :param a_matrix:
    :param v_matrix:
    :param find_model:
    :return:
    """
    # TODO: maybe give P, T and iterate here, so we simplify once (if it's a bottleneck)

    func_expression, symbols_to_variables_dict = None, None
    cur_func = G.vertices[i].function
    if cur_func is None:
        raise AssertionError("Nodes with unset functions shouldn't be handled here")
    if isinstance(cur_func, sympy.FunctionClass):
        # instantiate the function expression
        arg_symbols = [sympy.symbols("x_{}".format(j)) for j in range(len(predecessors_vars[i, p, t]))]
        func_expression = cur_func(*arg_symbols)
    elif isinstance(cur_func, logic.BooleanSymbolicFunc):
        func_expression = cur_func.formula
        arg_symbols = cur_func.input_vars
    else:
        # try a constant function
        try:
            if cur_func(None) in [True, False, sympy.true, sympy.false]:
                func_expression = sympy.true if cur_func(None) else sympy.false
                arg_symbols = [None] * len(G.vertices[i].predecessors())
            else:
                raise ValueError("Unkown type of function - " + str(type(cur_func)))
        except TypeError:
            raise ValueError("Unkown type of function (non-constant lambda functions aren't allowed)")
    simplified_formula = sympy.simplify(func_expression)
    symbols_to_variables_dict = dict(zip(arg_symbols, predecessors_vars[i, p, t]))
    func_var = build_logic_function_vars(simplified_formula, model, "Consistency_{}_".format(i),
                                         symbols_to_variables_dict)
    model.addConstr(v_matrix[i, p, t + 1] >= func_var + a_matrix[p, t] - 1,
                    name="consistent_>=_{}_{}_{}".format(i, p, t))
    model.addConstr(v_matrix[i, p, t + 1] <= func_var - a_matrix[p, t] + 1,
                    name="consistent_<=_{}_{}_{}".format(i, p, t))


# noinspection PyArgumentList
def direct_graph_to_ilp_with_keys(G, max_len=None, max_num=None,
                                  model_type_restriction=FunctionTypeRestriction.NONE):
    total_start = time.time()
    part_start = time.time()
    T = 2**len(G.vertices) if not max_len else max_len
    P = 2**len(G.vertices) if not max_num else max_num
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
    state_keys = [[unique_state_keys([v_matrix[i, p, t] for i in range(n)]) for t in range(T+1)]
                  for p in range(P)]
    predecessors_vars = numpy.array([[[[v_matrix[vertex.index, p, t] for vertex in G.vertices[i].predecessors()]
                                     for t in range(T+1)] for p in range(P)] for i in range(n)])
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
                for p, t in itertools.product(range(P), range(T)):
                    if simplify_general_boolean and not find_model:
                        add_simplified_consistency_constraints(G, model, i, p, t, predecessors_vars,
                                                               a_matrix, v_matrix)
                    else:
                        add_truth_table_consistency_constraints(G, model, i, p, t, predecessors_vars,
                                                                a_matrix, v_matrix, find_model_f_vars)
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
                                                               upper_bound=2**unique_key_slicing_size,
                                                               name_prefix="simple>_{}_{}".format(p, t))
        model.addConstr(strictly_larger_ind >= a_matrix[p, t],
                        name="simple_{}_{}".format(p, t))
    # model.update()

    # print "Time taken for simple constraints preparation:{:.2f} seconds".format(time.time() - part_start)
    # part_start = time.time()

    # UNIQUE
    # To reduce symmetry, using the order defined by a unique state id function,
    # constraint each attractor to have its final state be the largest one, and
    # constraint the final states of attractors to be monotone decreasing amongst.
    # Also requires the k active attractors to be the first ones. SHOULD TODO: verify
    # result in only one optimal solution.
    for p in range(P-1):
        # as long as keys are uniques, this forces uniqueness if p and p + 1 are both active
        strictly_larger_ind = create_state_keys_comparison_var(model=model,
                                                               first_state_keys=state_keys[p + 1][T-1],
                                                               second_state_keys=state_keys[p][T-1],
                                                               include_equality=False,
                                                               upper_bound=2**unique_key_slicing_size,
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
    return model, [a_matrix[p, T-1] for p in range(P)]


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
    n_attracotrs = int(round(model.ObjVal))
    if n_attracotrs != model.ObjVal:
        print "warning, model solved with non-integral objective value {}".format(model.ObjVal)
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
    assert len([None for p in range(P) if int(round(a_variables[p][T-1].X)) == 1]) == n_attracotrs
    for attractor_number in range(int(n_attracotrs)):
        p = P - attractor_number - 1
        length = len([None for t in range(T) if int(round(a_variables[p][t].X)) == 1])
        print "Attractor #{}, length {}".format(attractor_number + 1, length)
        # TODO: support for original graph names?
        # print reduce(lambda a, b: "{}\t{}".format(a, b), ["v_{}".format(i) for i in range(n)])
        for t in range(T - length, T):
            # noinspection PyTypeChecker
            print reduce(lambda a, b: "{}{}".format(a, b), [int(round(v_variables[i][p][t].X)) for i in range(n)])


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
