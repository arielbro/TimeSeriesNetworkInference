import gurobipy
import sympy
import numpy
import itertools


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


def unique_state_key(ordered_state_variables):
    """
    For small graphs, assign a unique numbering to a state, by summing 2**i over active vertices.
    :param Ordered_state_variables: an ordered fixed iterable of vertex state variables.
    :return: A key identifying this state uniquely among all other states.
    """
    if len(ordered_state_variables) > 29:
        raise Exception("Can't use unique state key with graphs of size >=30")
    return sum(2**i * var for i, var in enumerate(ordered_state_variables))


def direct_graph_to_ilp(G, max_len=None, max_num=None, find_bool_model=False):
    T = 2**len(G.vertices) if not max_len else max_len
    P = 2**len(G.vertices) if not max_num else max_num

    model = gurobipy.Model()

    a_matrix = numpy.matrix([[model.addVar(vtype=gurobipy.GRB.BINARY, name="a_{}_{}".format(p, t))
                              for t in range(T+1)] for p in range(P)])
    v_matrix = numpy.array([[[model.addVar(vtype=gurobipy.GRB.BINARY, name="v_{}_{}_{}".format(i, p, t))
                               for t in range(T+1)] for p in range(P)] for i in range(len(G.vertices))])
    model.update()
    predecessors_vars = lambda i, p, t: [v_matrix[vertex.index, p, t] for vertex in G.vertices[i].predecessors()]

    if find_bool_model:
        boolean_vars_dict = dict()
        for i in range(len(G.vertices)):
            for j in range(len(G.vertices[i].predecessors())):
                new_var = model.addVar(vtype=gurobipy.GRB.BINARY, name="f_{}_{}".format(i, j))
                model.update()
                boolean_vars_dict[(i, j)] = new_var

    for i, p, t in itertools.product(range(len(G.vertices)), range(P), range(T + 1)):
        # assert activity var meaning (ACTIVITY_SWITCH, MONOTONE, IF_NON_ACTIVE)
        model.addConstr(a_matrix[p, t] >= v_matrix[i, p, t], name="activity_switch_{}_{}_{}".format(i, p, t))
        if t != T:
            model.addConstr(a_matrix[p, t + 1] >= a_matrix[p, t], name="monotone_{}_{}_{}".format(i, p, t))
        else:
            model.addConstr(a_matrix[p, T] <= a_matrix[p, T - 1], name="if_non_active_{}_{}".format(i, p))
        model.update()


        # assert consistent and stable
        in_degree = len(G.vertices[i].predecessors())
        if in_degree == 0 and t < T:  # stable, i.e. a[i,p,t] => v[i,p,t+1] = v[i, p , t].
            model.addConstr(v_matrix[i, p, t + 1] >= v_matrix[i, p, t] - (1 - a_matrix[p, t]),
                            name="stable_>=_{}_{}_{}".format(i, p, t))
            model.addConstr(v_matrix[i, p, t + 1] <= v_matrix[i, p, t] + (1 - a_matrix[p, t]),
                            name="stable_<=_{}_{}_{}".format(i, p, t))
            model.update()
        elif t < T:
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
                desired_val = boolean_vars_dict[i, var_comb_index] if find_bool_model else \
                    (1 if G.vertices[i].function(*var_combination) else 0)
                # a[p, t] & indicator => v[i,p,t+1] = f(var_combination).
                # For x&y => a=b, require a <= b + (2 -x -b), a >= b - (2 -x -y)
                model.addConstr(v_matrix[i, p, t + 1] >= desired_val - (2 - indicator - a_matrix[p, t]),
                                name="consistent_>=_{}_{}_{}".format(i, p, t))
                model.addConstr(v_matrix[i, p, t + 1] <= desired_val + (2 - indicator - a_matrix[p, t]),
                                name="consistent_<=_{}_{}_{}".format(i, p, t))
        model.update()
        # cyclic constraints
        # a[p, 0] >> v[i, p, 0] == v[i, p, T]
        if t == 0:
            model.addConstr(v_matrix[i, p, 0] >= v_matrix[i, p, T] - (1 - a_matrix[p, 0]),
                            name="cyclic_>=_{}_{}_{}".format(i, p, t))
            model.addConstr(v_matrix[i, p, 0] <= v_matrix[i, p, T] + (1 - a_matrix[p, 0]),
                            name="cyclic_<=_{}_{}_{}".format(i, p, t))
        elif t < T:  # (~a[p, t - 1] & a[p, t]) >> v[i, p, t] = v[i, p, T], assumes monotonicity of a_matrix
            model.addConstr(v_matrix[i, p, t] >= v_matrix[i, p, T] - (1 - a_matrix[p, t] + a_matrix[p, t - 1]),
                            name="cyclic_>=_{}_{}_{}".format(i, p, t))
            model.addConstr(v_matrix[i, p, t] <= v_matrix[i, p, T] + (1 - a_matrix[p, t] + a_matrix[p, t - 1]),
                            name="cyclic_<=_{}_{}_{}".format(i, p, t))

    # SIMPLE
    for t, p in itertools.product(range(1, T), range(P)):
        # (a[p, t] & a[p, t-1]) >> ~EQ(p, p, t, T)
        equality_indicator_vars = [model.addVar(vtype=gurobipy.GRB.BINARY,
                                   name="eq_simple_ind_{}_{}_{}".format(i, p, t)) for i in range(len(G.vertices))]
        for i in range(len(G.vertices)):
            model.addConstr(equality_indicator_vars[i] <= 1 + v_matrix[i, p, t] - v_matrix[i, p, T],
                            name="eq_simple_ind_0_{}_{}_{}".format(i, p, t))
            model.addConstr(equality_indicator_vars[i] <= 1 - v_matrix[i, p, t] + v_matrix[i, p, T],
                            name="eq_simple_ind_1_{}_{}_{}".format(i, p, t))
            model.addConstr(equality_indicator_vars[i] >= -1 + v_matrix[i, p, t] + v_matrix[i, p, T],
                            name="eq__simple_ind_2_{}_{}_{}".format(i, p, t))
            model.addConstr(equality_indicator_vars[i] >= 1 - v_matrix[i, p, t] - v_matrix[i, p, T],
                            name="eq_simple_ind_3_{}_{}_{}".format(i, p, t))
        # it holds that ~EQ(p, p, t, T) <=> sum(equality_indicator_vars) <  len(G.vertices), now create the >> part
        model.addConstr(sum(equality_indicator_vars) <= len(G.vertices) + 1 - a_matrix[p, t] - a_matrix[p, t - 1],
                        name="simple_{}_{}".format(t, p))

    # UNIQUE
    # for t, (p1, p2) in itertools.product(range(T), itertools.combinations(range(P), 2)):
    #     # (a[p1, T] & a[p2, t]) >> ~EQ(p1, p2, T, t)
    #     equality_indicator_vars = [model.addVar(vtype=gurobipy.GRB.BINARY,
    #                                name="eq_unique_ind_{}_{}_{}_{}".format(i, p1, p2, t))
    #                                for i in range(len(G.vertices))]
    #     for i in range(len(G.vertices)):
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

    # To reduce symmetry, using the order defined by a unique state id function,
    # constraint each attractor to have its final state be the largest one, and
    # constraint the final states of attractors to be monotone decreasing amongst.
    # Also requires the k active attractors to be the first ones. SHOULD TODO: verify
    # result in only one optimal solution.
    final_states_keys = [unique_state_key([v_matrix[i, p, T] for i in range(len(G.vertices))])
                         for p in range(P)]
    for p in range(P):
        for t in range(T):
            current_state_key = unique_state_key([v_matrix[i, p, t] for i in range(len(G.vertices))])
            # works for inactive states attractors/time points too!
            model.addConstr(final_states_keys[p] >= current_state_key,
                            name="key_order_in_attractor_{}_{}".format(p, t))

        if p != P - 1:  # as long as keys are uniques, this forces uniqueness if p and p + 1 are both active
            model.addConstr(final_states_keys[p] >= final_states_keys[p + 1]
                            - 1 + a_matrix[p, T] + a_matrix[p + 1, T],
                                                            name="key_order_between_attractors_{}".format(p))
    model.update()

    # Constraint the number of active attractors using 2**#input_nodes, P as lower and upper bounds.
    # lower bound can only be used if the maximal theoretical attractor length is allowed.
    n_inputs = len([v for v in G.vertices if len(v.predecessors()) == 0])
    model.addConstr(sum(a_matrix[p, T] for p in range(P)) <= P, name="upper_objective_bound")
    if T == 2**len(G.vertices):
        model.addConstr(sum(a_matrix[p, T] for p in range(P)) >= 2**n_inputs, name="lower_objective_bound")
    model.update()

    # print_model_constraints(model)
    # print model
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


def print_model_constraints(model):
    quadruples = list(get_matrix_coos(model))
    constr_attrs = [(constr.Sense, constr.RHS) for constr in model.getConstrs()]
    for constraint_index in range(max(row_index for row_index, _, _, _ in quadruples)):
        con_str = None
        for row_index, constr_name, var_name, coeff in quadruples:
            if row_index == constraint_index:
                if not con_str:
                    con_str = constr_name + ": "
                con_str += " {}{}{}".format(str(coeff) if coeff not in [1.0, -1.0] else "",
                                            "-" if coeff == -1.0 else "+" if coeff == 1 else "", var_name)
        con_str += " {} {}".format(*constr_attrs[constraint_index])
        print con_str


