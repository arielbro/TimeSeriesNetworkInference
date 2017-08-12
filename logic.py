import sympy
import itertools
import numpy

def pre_randomized_boolean_symbolic_func(boolean_outputs, *variables):
    formula = sympy.true
    for b_output, terms in zip(boolean_outputs, itertools.product(*[[~var, var] for var in variables])):
        if b_output: # TODO: Karnaugh maps?
            formula = formula | sympy.And(*terms)
    return formula


def pre_randomized_boolean_func(boolean_outputs):
    input_variables = [sympy.symbols("x" + str(i)) for i in range(len(boolean_outputs))]
    symbolic_func = pre_randomized_boolean_symbolic_func(boolean_outputs, *input_variables)
    return lambda *boolean_inputs: \
        symbolic_func.subs({input_variables[i]: boolean_inputs[i] for i in range(len(boolean_inputs))})


# TODO: create a pretty print option (e.g. a callable class creating partials with a homemade str()
def variable_bound_boolean_func(function_output_vars, *variables):
    formula = sympy.false
    for output_var, terms in zip(function_output_vars, itertools.product(*[[~var, var] for var in variables])):
        formula = formula | (sympy.And(*terms) & output_var)
    return formula


def formula_length(formula):
    # defined as the number of (non-unique) atoms in the formula
    if formula.is_Atom:
        return 1
    else:
        return sum(formula_length(arg) for arg in formula.args)


def get_attractors_formula(G, P, T):
    a_matrix = numpy.matrix([[sympy.symbols("a_{}_{}".format(p, t)) for t in range(T+1)] for p in range(P)])
    v_matrix = numpy.array([[[sympy.symbols("v_{}_{}_{}".format(i, p, t)) for t in range(T+1)] for p in range(P)]
                             for i in range(len(G.vertices))])

    ACTIVITY_SWITCH = lambda p: sympy.And(*[~a_matrix[p, t] >> sympy.And(*[~v_matrix[i, p, t]
                                            for i in range(len(G.vertices))]) for t in range(T + 1)])
    MONOTONE = lambda p: sympy.And(*[a_matrix[p, t] >> a_matrix[p, t+1] for t in range(T)])
    IF_NON_ACTIVE = lambda p: ~a_matrix[p, T-1] >> ~a_matrix[p, T]
    # ACTIVE = lambda p: a_matrix[p, T-1]


    predecessors_vars = lambda i, p, t: [v_matrix[vertex.index, p, t] for vertex in G.vertices[i].predecessors()]

    CONSISTENT = lambda p: sympy.And(*[sympy.And(*[
                                     a_matrix[p, t] >> (sympy.Equivalent(v_matrix[i, p, t+1],
                                     G.vertices[i].function(*predecessors_vars(i, p, t))))
                                     for i in range(len(G.vertices)) if len(G.vertices[i].predecessors()) > 0])
                                     for t in range(T)])

    STABLE = lambda p: sympy.And(*[sympy.And(*[a_matrix[p, t] >>
                                               sympy.Equivalent(v_matrix[i, p, t], v_matrix[i, p, t+1])
                                               for t in range(T)]) for i in range(len(G.vertices)) if
                                               len(G.vertices[i].predecessors()) == 0])

    EQ = lambda p1, p2, t1, t2: sympy.And(*[sympy.Equivalent(v_matrix[i, p1, t1], v_matrix[i, p2, t2])
                                            for i in range(len(G.vertices))])
    CYCLIC = lambda p: (a_matrix[p, 0] >> EQ(p, p, 0, T)) & \
                       (sympy.And(*[(~a_matrix[p, t - 1] & a_matrix[p, t]) >> EQ(p, p, t, T)
                                                                       for t in range(1, T + 1)]))
    SIMPLE = lambda p: sympy.And(*[(a_matrix[p, t] & a_matrix[p, t-1]) >> ~EQ(p, p, t, T) for t in range(1, T)])

    UNIQUE = lambda p1: sympy.And(*[sympy.And(*[(a_matrix[p1, T] & a_matrix[p2, t]) >> ~EQ(p1, p2, T, t)
                                                for p2 in range(p1 + 1, P)]) for t in range(T)])

    ATTRACTORS = sympy.And(*[ACTIVITY_SWITCH(p) & MONOTONE(p) & IF_NON_ACTIVE(p) & CONSISTENT(p) &
                             STABLE(p) & CYCLIC(p) & SIMPLE(p)
                             & UNIQUE(p) for p in range(P)])

    # print ACTIVITY_SWITCH(0)
    # print MONOTONE(0)
    # print IF_NON_ACTIVE(0)
    # print CONSISTENT(0)
    # print STABLE(0)
    # print CYCLIC(0)
    # print SIMPLE(0)
    # print UNIQUE(0)

    return ATTRACTORS, [a_matrix[p, T] for p in range(P)]  #, a_matrix, v_matrix


def get_attractorlb_lengthub_formula(G, P, T):
    ATTRACTORS, activity_formulas = get_attractors_formula(G, P, T)
    ATTRACTORS = sympy.And(*([ATTRACTORS] + activity_formulas))
    return ATTRACTORS  #, a_matrix, v_matrix
