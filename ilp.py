import gurobipy
import sympy
import numpy
import itertools
import logic
import time
import math
from graphs import FunctionTypeRestriction

unique_key_slicing_size = 29  # TODO: find elegant reformatting for this


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


def unique_state_keys(model, ordered_state_variables, name_prefix):
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
        model.addConstr((M + 1)*z >= a - b + last_var, name="{}_{}_>=constraint".format(name_prefix, i))
        model.addConstr((M + 1)*z <= a - b + last_var + M, name="{}_{}_<=constraint".format(name_prefix, i))
        model.update()
        last_var = z
        # print "a_{}={}, b_{}={}, M={}, M+1={}".format(len(first_state_keys) - i-1, a, len(first_state_keys) -i-1, b, M, M+1)
    return last_var


# noinspection PyArgumentList
def direct_graph_to_ilp(G, max_len=None, max_num=None, find_model=False,
                        model_type_restriction=FunctionTypeRestriction.NONE):
    start = time.time()
    T = 2**len(G.vertices) if not max_len else max_len
    P = 2**len(G.vertices) if not max_num else max_num
    n = len(G.vertices)

    model = gurobipy.Model()

    a_matrix = numpy.matrix([[model.addVar(vtype=gurobipy.GRB.BINARY, name="a_{}_{}".format(p, t))
                              for t in range(T+1)] for p in range(P)]) # TODO: fill for t=-1, and simplify code
    v_matrix = numpy.array([[[model.addVar(vtype=gurobipy.GRB.BINARY, name="v_{}_{}_{}".format(i, p, t))
                               for t in range(T+1)] for p in range(P)] for i in range(n)])
    model.update()

    state_keys = [[unique_state_keys(model=model, ordered_state_variables=[v_matrix[i, p, t] for i in range(n)],
                                     name_prefix="state_keys_{}_{}".format(p, t)) for t in range(T+1)]
                  for p in range(P)]

    predecessors_vars = lambda i, p, t: [v_matrix[vertex.index, p, t] for vertex in G.vertices[i].predecessors()]
    if find_model:
        if model_type_restriction == FunctionTypeRestriction.NONE:
            truth_table_vars_dict = dict()
            for i in range(n):
                for truth_table_index in range(2**len(G.vertices[i].predecessors())):
                    new_var = model.addVar(vtype=gurobipy.GRB.BINARY,
                                           name="f_{}_{}".format(i, truth_table_index))
                    model.update()
                    truth_table_vars_dict[(i, truth_table_index)] = new_var
        elif model_type_restriction == FunctionTypeRestriction.SYMMETRIC_THRESHOLD\
                or model_type_restriction == FunctionTypeRestriction.SIMPLE_GATES:
            signs_threshold_vars_dict = dict()
            for i in range(n):
                n_inputs = len(G.vertices[i].predecessors())
                signs = []
                for input in range(n_inputs):
                    sign_var = model.addVar(vtype=gurobipy.GRB.BINARY, name="f_{}_signs_{}".format(i, input))
                    signs.append(sign_var)
                if model_type_restriction == FunctionTypeRestriction.SYMMETRIC_THRESHOLD:
                    threshold_expression = model.addVar(vtype=gurobipy.GRB.INTEGER, lb=0, ub=n_inputs+1,
                                             name="f_{}_threshold".format(i))
                else:
                    threshold_indicator = model.addVar(
                        vtype=gurobipy.GRB.BINARY, name="f_{}_gate_indicator".format(i))
                    threshold_expression = 1 + threshold_indicator * (n_inputs - 1)
                model.update()
                signs_threshold_vars_dict[i] = signs, threshold_expression

    for p, t in itertools.product(range(P), range(T)):
        # assert activity var meaning (ACTIVITY_SWITCH, MONOTONE, IF_NON_ACTIVE)
        # for i in range(n):
        #     model.addConstr(a_matrix[p, t] >= v_matrix[i, p, t], name="activity_switch_{}_{}_{}".format(i, p, t))
        model.addConstr(n * a_matrix[p, t] >= sum(v_matrix[i, p, t] for i in range(n)),
                        name="activity_{}_{}".format(p, t))
        model.addConstr(a_matrix[p, t + 1] >= a_matrix[p, t], name="monotone_{}_{}".format(p, t))
        model.update()

    for p in range(P):
        model.addConstr(a_matrix[p, T] <= a_matrix[p, T - 1], name="if_non_active_{}".format(p))
        if p != P - 1:
            model.addConstr(a_matrix[p, T] <= a_matrix[p + 1, T], name="active_order_{}".format(p))

    for i, p, t in itertools.product(range(n), range(P), range(T + 1)):

        # assert consistent and stable
        in_degree = len(G.vertices[i].predecessors())
        if in_degree == 0 and t < T:  # stable, i.e. a[i,p,t] => v[i,p,t+1] = v[i, p , t].
            model.addConstr(v_matrix[i, p, t + 1] >= v_matrix[i, p, t] - (1 - a_matrix[p, t]),
                            name="stable_>=_{}_{}_{}".format(i, p, t))
            model.addConstr(v_matrix[i, p, t + 1] <= v_matrix[i, p, t] + (1 - a_matrix[p, t]),
                            name="stable_<=_{}_{}_{}".format(i, p, t))
            model.update()
        elif t < T:
            if isinstance(G.vertices[i].function, logic.SymmetricThresholdFunction) or \
                    (find_model and model_type_restriction != FunctionTypeRestriction.NONE):
                if find_model and model_type_restriction != FunctionTypeRestriction.NONE:
                    signs, threshold = signs_threshold_vars_dict[i]
                else:
                    signs = G.vertices[i].function.signs
                    threshold = G.vertices[i].function.threshold
                input_sum_expression = 0
                for predecessors_index, (sign, var) in enumerate(zip(signs, predecessors_vars(i, p, t))):
                    if model_type_restriction != FunctionTypeRestriction.NONE:
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
                        model.update()
                        input_sum_expression += z
                    else:
                        input_sum_expression += var if sign else 1 - var  # sign is bool, unintuitive?
                # input_sum_expression >= threshold <=> v_matrix[i, p, t + 1] = 1, if activity is on.
                model.addConstr((in_degree + 1)*(v_matrix[i, p, t + 1] + 1 - a_matrix[p, t]) >=
                                (input_sum_expression - threshold + 1),
                                name="monotone_func_consistency_>=_{}_{}_{}".format(i, p, t))
                model.addConstr((in_degree + 1)*(v_matrix[i, p, t + 1] - 1 + a_matrix[p, t]) <=
                                (in_degree + input_sum_expression - threshold + 1),
                                name="monotone_func_consistency_<=_{}_{}_{}".format(i, p, t))
                model.update()
            else:
                for var_comb_index, var_combination in enumerate(itertools.product((False, True), repeat=in_degree)):
                    # this expression is |in_degree| iff their states agrees with var_combination
                    indicator_expression = sum(v if state else 1 - v for (v, state) in
                                               zip(predecessors_vars(i, p, t), var_combination))
                    # TODO: see if the indicator expression can be used without the matching variable
                    indicator = model.addVar(vtype=gurobipy.GRB.BINARY,
                                             name="truth_table_row_indicator_{}_{}_{}_{}".format(i, p, t, var_comb_index))
                    model.addConstr(in_degree * indicator <= indicator_expression,
                                    name="truth_table_row_indicator_<=_{}_{}_{}_{}".format(i, p, t, var_comb_index))
                    model.addConstr(indicator >= indicator_expression - in_degree + 1,
                                    name="truth_table_row_indicator_>=_{}_{}_{}_{}".format(i, p, t, var_comb_index))
                    # noinspection PyUnboundLocalVariable
                    desired_val = truth_table_vars_dict[i, var_comb_index] if find_model else \
                        (1 if G.vertices[i].function(*var_combination) else 0)
                    # a[p, t] & indicator => v[i,p,t+1] = f(var_combination).
                    # For x&y => a=b, require a <= b + (2 -x -b), a >= b - (2 -x -y)
                    model.addConstr(v_matrix[i, p, t + 1] >= desired_val - (2 - indicator - a_matrix[p, t]),
                                    name="consistent_>=_{}_{}_{}".format(i, p, t))
                    model.addConstr(v_matrix[i, p, t + 1] <= desired_val + (2 - indicator - a_matrix[p, t]),
                                    name="consistent_<=_{}_{}_{}".format(i, p, t))
                    model.update()

    # CYCLIC
    for p, t in itertools.product(range(P), range(T)):
        # (~a[p, t - 1] & a[p, t]) >> EQ(p, p, t, T), assumes monotonicity of a_matrix (a[p, -1] assumed 0)
        larger_ind = create_state_keys_comparison_var(model=model,
                                                      first_state_keys=state_keys[p][T],
                                                      second_state_keys=state_keys[p][t],
                                                      include_equality=True,
                                                      upper_bound=2**unique_key_slicing_size,
                                                      name_prefix="cyclic>_{}_{}".format(p, t).
                                                      format(p))
        smaller_ind = create_state_keys_comparison_var(model=model,
                                                       first_state_keys=state_keys[p][t],
                                                       second_state_keys=state_keys[p][T],
                                                       include_equality=True,
                                                       upper_bound=2**unique_key_slicing_size,
                                                       name_prefix="cyclic<_{}_{}".format(p, t).
                                                       format(p))
        last_activity = a_matrix[p, t - 1] if t > 0 else 0
        model.addConstr(larger_ind + smaller_ind >= 1 - last_activity + a_matrix[p, t],
                        name="cyclic_{}_{}".format(p, t))
        # for i in range(n):
        #     model.addConstr(v_matrix[i, p, t] >= v_matrix[i, p, T] + a_matrix[p, t] - last_activity - 1,
        #                     name="cyclic<_{}_{}_{}".format(i, p, t))
        #     model.addConstr(v_matrix[i, p, t] <= v_matrix[i, p, T] - a_matrix[p, t] + last_activity + 1,
        #                     name="cyclic>_{}_{}_{}".format(i, p, t))
    model.update()

    # SIMPLE
    # for t, p in itertools.product(range(1, T), range(P)):
    #     # (a[p, t] & a[p, t-1]) >> ~EQ(p, p, t, T)
    #     equality_indicator_vars = [model.addVar(vtype=gurobipy.GRB.BINARY,
    #                                name="eq_simple_ind_{}_{}_{}".format(i, p, t)) for i in range(n)]
    #     for i in range(n):
    #         model.addConstr(equality_indicator_vars[i] <= 1 + v_matrix[i, p, t] - v_matrix[i, p, T],
    #                         name="eq_simple_ind_0_{}_{}_{}".format(i, p, t))
    #         model.addConstr(equality_indicator_vars[i] <= 1 - v_matrix[i, p, t] + v_matrix[i, p, T],
    #                         name="eq_simple_ind_1_{}_{}_{}".format(i, p, t))
    #         model.addConstr(equality_indicator_vars[i] >= -1 + v_matrix[i, p, t] + v_matrix[i, p, T],
    #                         name="eq__simple_ind_2_{}_{}_{}".format(i, p, t))
    #         model.addConstr(equality_indicator_vars[i] >= 1 - v_matrix[i, p, t] - v_matrix[i, p, T],
    #                         name="eq_simple_ind_3_{}_{}_{}".format(i, p, t))
    #     # it holds that ~EQ(p, p, t, T) <=> sum(equality_indicator_vars) <  len(G.vertices), now create the >> part
    #     model.addConstr(sum(equality_indicator_vars) <= len(G.vertices) + 1 - a_matrix[p, t] - a_matrix[p, t - 1],
    #                     name="simple_{}_{}".format(t, p))
    # model.update()

    for t, p in itertools.product(range(1, T), range(P)):
        # (a[p, t] & a[p, t-1]) >> ~EQ(p, p, t, T)
        strictly_larger_ind = create_state_keys_comparison_var(model=model, first_state_keys=state_keys[p][T],
                                                               second_state_keys=state_keys[p][t],
                                                               include_equality=False,
                                                  upper_bound= 2**unique_key_slicing_size,
                                                               name_prefix="simple>_{}_{}".format(p, t))
        strictly_smaller_ind = create_state_keys_comparison_var(model=model, first_state_keys=state_keys[p][t],
                                                                second_state_keys=state_keys[p][T],
                                                                include_equality=False,
                                                  upper_bound= 2**unique_key_slicing_size,
                                                                name_prefix="simple<_{}_{}".format(p, t))
        model.addConstr(strictly_larger_ind + strictly_smaller_ind - a_matrix[p, t] - a_matrix[p, t-1] >= -1,
                        name="simple_{}_{}".format(p, t))
    model.update()

    # UNIQUE
    # for t, (p1, p2) in itertools.product(range(T), itertools.combinations(range(P), 2)):
    #     # (a[p1, T] & a[p2, t]) >> ~EQ(p1, p2, T, t)
    #     equality_indicator_vars = [model.addVar(vtype=gurobipy.GRB.BINARY,
    #                                name="eq_unique_ind_{}_{}_{}_{}".format(i, p1, p2, t))
    #                                for i in range(n)]
    #     for i in range(n):
    #         model.addConstr(equality_indicator_vars[i] <= 1 + v_matrix[i, p1, T] - v_matrix[i, p2, t],
    #                         name="eq_unique_ind_0_{}_{}_{}_{}".format(i, p1, p2, t))
    #         model.addConstr(equality_indicator_vars[i] <= 1 - v_matrix[i, p1, T] + v_matrix[i, p2, t],
    #                         name="eq_unique_ind_1_{}_{}_{}_{}".format(i, p1, p2, t))
    #         model.addConstr(equality_indicator_vars[i] >= -1 + v_matrix[i, p1, T] + v_matrix[i, p2, t],
    #                         name="eq_unique_ind_2_{}_{}_{}_{}".format(i, p1, p2, t))
    #         model.addConstr(equality_indicator_vars[i] >= 1 - v_matrix[i, p1, T] - v_matrix[i, p2, t],
    #                         name="eq_unique_ind_3_{}_{}_{}_{}".format(i, p1, p2, t))
    #     # it holds that ~EQ(p1, p2, T, t) <=> sum(equality_indicator_vars) <  len(G.vertices), now create the >> part
    #     model.addConstr(sum(equality_indicator_vars) <= len(G.vertices) + 1 - a_matrix[p1, T] - a_matrix[p2, t],
    #                     name="unique_{}_{}_{}".format(p1, p2, t))
    # model.update()

    # To reduce symmetry, using the order defined by a unique state id function,
    # constraint each attractor to have its final state be the largest one, and
    # constraint the final states of attractors to be monotone decreasing amongst.
    # Also requires the k active attractors to be the first ones. SHOULD TODO: verify
    # result in only one optimal solution.
    for p in range(P):
        for t in range(T):
            larger_ind = create_state_keys_comparison_var(model=model, first_state_keys=state_keys[p][T],
                                                          second_state_keys=state_keys[p][t],
                                                          include_equality=True,
                                                          upper_bound=2**unique_key_slicing_size,
                                                          name_prefix="key_order_in_attractor_{}_{}".format(p, t))
            # works for inactive states attractors/time points too!
            model.addConstr(larger_ind >= 1, name="key_order_in_attractor_{}_{}".format(p, t))

        if p != P - 1:  # as long as keys are uniques, this forces uniqueness if p and p + 1 are both active
            strictly_larger_ind = create_state_keys_comparison_var(model=model,
                                                                   first_state_keys=state_keys[p + 1][T],
                                                                   second_state_keys=state_keys[p][T],
                                                                   include_equality=False,
                                                                   upper_bound=2**unique_key_slicing_size,
                                                                   name_prefix="key_order_between_attractors_{}".
                                                                   format(p))
            model.addConstr(strictly_larger_ind - a_matrix[p, T] - a_matrix[p + 1, T] >= -1,
                            name="key_order_between_attractors_{}".format(p))
        model.update()

    # Constraint the number of active attractors using 2**#input_nodes, P as lower and upper bounds.
    # lower bound can only be used if the maximal theoretical attractor length is allowed.
    n_inputs = len([v for v in G.vertices if len(v.predecessors()) == 0])
    model.addConstr(sum(a_matrix[p, T] for p in range(P)) <= P, name="upper_objective_bound")
    if T >= 2**len(G.vertices):
        model.addConstr(sum(a_matrix[p, T] for p in range(P)) >= min(2**n_inputs, P), name="lower_objective_bound")
    model.update()

    # print_model_constraints(model)
    # print model
    print "Time taken for model preparation:{:.2f} seconds".format(time.time() - start)
    return model, [a_matrix[p, T] for p in range(P)]


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
        length = len([None for t in range(T) if a_variables[p][t].X == 1])
        print "Attractor #{}, length {}".format(attractor_number + 1, length)
        # TODO: support for original graph names?
        print reduce(lambda a, b: "{}\t{}".format(a, b), ["v_{}".format(i) for i in range(n)])
        for t in range(T - length, T):
            # noinspection PyTypeChecker
            print reduce(lambda a, b: "{}\t{}".format(a, b), [int(round(v_variables[i][p][t].X)) for i in range(n)])
        print "\n"


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
