import sympy
import itertools
import math
import numpy
from utility import list_repr


class BooleanSymbolicFunc:
    def __init__(self, boolean_outputs):
        self.truth_table_outputs = boolean_outputs
        # assumes boolean_outputs is a power of 2
        n_inputs = math.frexp(len(boolean_outputs))[1] - 1
        boolean_inputs = [sympy.symbols("x_{}".format(i)) for i in range(n_inputs)]
        self.input_vars = boolean_inputs
        if n_inputs == 0:
            assert len(boolean_outputs) == 1
            self.formula = boolean_outputs[0]
            return
        formula = sympy.false
        for b_output, terms in zip(boolean_outputs, itertools.product(*[[~var, var] for var in boolean_inputs])):
            if b_output:  # TODO: Karnaugh maps?
                formula = formula | sympy.And(*terms)
        self.formula = formula

    def __call__(self, *input_values):
        if isinstance(self.formula, bool) or len(self.formula.free_symbols) == 0:
            return self.formula
        return self.formula.subs(zip(self.input_vars, input_values))

    def __str__(self):
        return " " + str(self.formula)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        try:
            for input_comb in itertools.product([False, True], repeat=len(self.input_vars)):
                if self(*input_comb) != other(*input_comb):
                    return False
        except ValueError:
            return False
        return True

    def __hash__(self):
        output_string = ""
        for input_comb in itertools.product([False, True], repeat=len(self.input_vars)):
            output_string += "1" if self(*input_comb) == True else "0"
        return hash(output_string)

    def __ne__(self, other):
        return not self == other


class SymmetricThresholdFunction:
    # TODO: implement in ILP model finding (threshold is not boolean, not supported there ATM)
    def __init__(self, signs, threshold):
        # translate signs to bool values, if not already
        self.signs = []
        for sign in signs:
            if isinstance(sign, bool):
                self.signs.append(sign)
            elif isinstance(sign, int):
                assert sign in [1, -1]
                self.signs.append(True if sign == 1 else False)
            else:
                raise ValueError("illegal type for signs:{}".format(type(sign)))
        self.threshold = threshold

    def __call__(self, *input_values):
        count = sum(1 if ((sign and val) or (not sign and not val)) else 0
                    for (sign, val) in zip(self.signs, input_values))
        return count >= self.threshold

    def __str__(self):
        return "signs={}, threshold={}".format(list_repr([1 if sign else -1 for sign in self.signs]),
                                               self.threshold)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        try:
            for input_comb in itertools.product([False, True], repeat=len(self.signs)):
                if self(*input_comb) != other(*input_comb):
                    return False
        except ValueError:
            return False
        return True

    def __hash__(self):
        output_string = ""
        for input_comb in itertools.product([False, True], repeat=len(self.signs)):
            output_string += "1" if self(*input_comb) == True else "0"
        return hash(output_string)

    def __ne__(self, other):
        return not self == other

    # TODO: optimize?
    @staticmethod
    def from_function(function, n_args):  # TODO: tests!!
        if n_args == 0:
            return None
        input_combinations = itertools.product([False, True], repeat=n_args)
        f_in_out_touples = {tuple(combination): function(*combination) for combination in input_combinations}
        signs = []
        # find signs
        for i in range(n_args):
            negative = False
            positive = False
            for combination in f_in_out_touples.keys():
                if not combination[i]:  # only need to check half
                    f_1 = f_in_out_touples[combination]
                    f_2 = f_in_out_touples[combination[:i] + (True,) + combination[i + 1:]]
                    if f_2 and not f_1:
                        positive = True
                    if f_1 and not f_2:
                        negative = True
            if positive and negative:
                raise ValueError("Tried to convert a non symmetric-threshold function")
            if not positive and not negative:
                # constant function
                assert len(set(f_in_out_touples.values())) == 1
                if True in set(f_in_out_touples.values()):
                    return SymmetricThresholdFunction(signs=[True] * n_args, threshold=0)
                else:
                    assert False in set(f_in_out_touples.values())
                    return SymmetricThresholdFunction(signs=[True] * n_args, threshold=n_args+1)
            else:
                signs.append(True if positive else False)

        # find out threshold
        threshold = None
        for i in range(1, n_args + 1):
            i_combs = [combination for combination in f_in_out_touples.keys() if
                       sum(1 for sign, val in zip(signs, combination) if (val == int(sign))) == i]
            outputs = set([f_in_out_touples[i_comb] for i_comb in i_combs])
            if len(outputs) != 1:
                raise ValueError("Tried to convert a non symmetric-threshold function")
            if outputs == {True}:
                if threshold is None:
                    threshold = i
            else:
                if threshold is not None:
                    raise ValueError("Tried to convert a non symmetric-threshold function")
        if threshold is None:
            raise ValueError("Tried to convert a non symmetric-threshold function")
        return SymmetricThresholdFunction(signs=signs, threshold=threshold)


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
                                                                       for t in range(1, T)]))
    SIMPLE = lambda p: sympy.And(*[(a_matrix[p, t] & a_matrix[p, t-1]) >> ~EQ(p, p, t, T) for t in range(1, T)])

    UNIQUE = lambda p1: sympy.And(*[sympy.And(*[(a_matrix[p1, T] & a_matrix[p2, t]) >> ~EQ(p1, p2, T, t)
                                                for p2 in range(p1 + 1, P)]) for t in range(T)])

    # to reduce symmetry
    ACTIVES_FIRST = lambda p: True if p == P - 1 else (~a_matrix[p, T] >> ~a_matrix[p + 1, T])

    ATTRACTORS = sympy.And(*[ACTIVITY_SWITCH(p) & MONOTONE(p) & IF_NON_ACTIVE(p) & CONSISTENT(p) &
                             STABLE(p) & CYCLIC(p) & SIMPLE(p)
                             & UNIQUE(p) & ACTIVES_FIRST(p) for p in range(P)])


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