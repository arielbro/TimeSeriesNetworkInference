from utility import list_repr
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


# def direct_graph_to_ilp(G, max_len=None, max_num=None, find_bool_model=False):
#     T = 2**len(G.vertices) if not max_len else max_len
#     P = 2**len(G.vertices) if not find_bool_model else find_bool_model
#
#     model = gurobipy.Model()
#
#     a_matrix = numpy.matrix([[model.addVar(vtype=gurobipy.GRB.BINARY, name="a_{}_{}".format(p, t))
#                               for t in range(T+1)] for p in range(P)])
#     v_matrix = numpy.matrix([[[model.addVar(vtype=gurobipy.GRB.BINARY, name="a_{}_{}_{}".format(i, p, t))
#                               for t in range(T+1)] for p in range(P)] for i in range(len(G.vertices))])
#
#     for i, p, t in itertools.product(range(len(G.vertices)), range(P), range(T + 1)):
#         # assert activity var meaning (ACTIVITY_SWITCH, MONOTONE, IF_NON_ACTIVE)
#         model.add a_matrix[p, t] = v_matrix[i, p, t]
#         ACTIVITY_SWITCH = lambda p: sympy.And(*[~a_matrix[p, t] >> sympy.And(*[~v_matrix[i, p, t]
#                                                 for i in range(len(G.vertices))]) for t in range(T + 1)])
#         MONOTONE = lambda p: sympy.And(*[a_matrix[p, t] >> a_matrix[p, t+1] for t in range(T)])
#         IF_NON_ACTIVE = lambda p: ~a_matrix[p, T-1] >> ~a_matrix[p, T]
#         # ACTIVE = lambda p: a_matrix[p, T-1]
#
#
#         predecessors_vars = lambda i, p, t: [v_matrix[vertex.index, p, t] for vertex in G.vertices[i].predecessors()]
#
#         CONSISTENT = lambda p: sympy.And(*[sympy.And(*[
#                                          a_matrix[p, t] >> (sympy.Equivalent(v_matrix[i, p, t+1],
#                                          G.vertices[i].function(*predecessors_vars(i, p, t))))
#                                          for i in range(len(G.vertices)) if len(G.vertices[i].predecessors()) > 0])
#                                          for t in range(T)])
#
#         STABLE = lambda p: sympy.And(*[sympy.And(*[a_matrix[p, t] >>
#                                                    sympy.Equivalent(v_matrix[i, p, t], v_matrix[i, p, t+1])
#                                                    for t in range(T)]) for i in range(len(G.vertices)) if
#                                                    len(G.vertices[i].predecessors()) == 0])
#
#         EQ = lambda p1, p2, t1, t2: sympy.And(*[sympy.Equivalent(v_matrix[i, p1, t1], v_matrix[i, p2, t2])
#                                                 for i in range(len(G.vertices))])
#         CYCLIC = lambda p: (a_matrix[p, 0] >> EQ(p, p, 0, T)) & \
#                            (sympy.And(*[(~a_matrix[p, t - 1] & a_matrix[p, t]) >> EQ(p, p, t, T)
#                                                                            for t in range(1, T + 1)]))
#         SIMPLE = lambda p: sympy.And(*[(a_matrix[p, t] & a_matrix[p, t-1]) >> ~EQ(p, p, t, T) for t in range(1, T)])
#
#         UNIQUE = lambda p1: sympy.And(*[sympy.And(*[(a_matrix[p1, T] & a_matrix[p2, t]) >> ~EQ(p1, p2, T, t)
#                                                     for p2 in range(p1 + 1, P)]) for t in range(T)])
#
#         ATTRACTORS = sympy.And(*[ACTIVITY_SWITCH(p) & MONOTONE(p) & IF_NON_ACTIVE(p) & CONSISTENT(p) &
#                                  STABLE(p) & CYCLIC(p) & SIMPLE(p)
#                                  & UNIQUE(p) for p in range(P)])
#
#     # print ACTIVITY_SWITCH(0)
#     # print MONOTONE(0)
#     # print IF_NON_ACTIVE(0)
#     # print CONSISTENT(0)
#     # print STABLE(0)
#     # print CYCLIC(0)
#     # print SIMPLE(0)
#     # print UNIQUE(0)

