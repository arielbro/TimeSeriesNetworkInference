import sympy
import itertools
import math
import numpy
from .utility import list_repr
import time
import random


class BooleanSymbolicFunc(object):
    def __init__(self, input_names=None, boolean_outputs=None, formula=None, simplify_boolean_outputs=False):
        # make all fields immutable, so the function can be shallow copied safely.
        self._boolean_outputs = None if boolean_outputs is None else tuple(boolean_outputs)

        if formula is not None:
            self.input_vars = tuple(sorted(formula.free_symbols, key=lambda x: x.name))
            self.formula = formula
            return

        if input_names is None:
            n_inputs = int(math.log(len(boolean_outputs), 2))
            input_names = ['input_{}'.format(i) for i in range(n_inputs)]

        if len(input_names) != math.frexp(len(boolean_outputs))[1] - 1:
            raise ValueError("non-matching length for variable names list and boolean outputs list")
        # self.truth_table_outputs = boolean_outputs
        # assumes boolean_outputs is a power of 2
        n_inputs = len(input_names)
        boolean_inputs = tuple(sympy.symbols(name) for name in input_names)
        self.input_vars = boolean_inputs
        if n_inputs == 0:
            assert len(boolean_outputs) == 1
            self.formula = boolean_outputs[0]
            return
        # TODO: Karnaugh maps? Sympy simplification?
        positive_row_clauses = [sympy.And(*terms) for b_output, terms in zip(
            boolean_outputs, itertools.product(*[[~var, var] for var in boolean_inputs])) if b_output]
        self.formula = sympy.Or(*positive_row_clauses)
        if simplify_boolean_outputs:
            start = time.time()
            self.formula = sympy.simplify(self.formula)

    @property
    def boolean_outputs(self):
        if self._boolean_outputs is None:
            self._boolean_outputs = tuple(self(*row) for row in itertools.product([False, True],
                                                                                  repeat=len(self.input_vars)))
        return self._boolean_outputs

    @property
    def formula(self):
        return self._formula

    @formula.setter
    def formula(self, value):
        self._boolean_outputs = None
        self._formula = value
        # print "set formula value: {}".format(value)
        # print "formula type: {}".format(type(value))
        if isinstance(value, sympy.Basic):
            # print "set lambdified formula"
            self._lambdified_formula = sympy.lambdify(self.input_vars, self.formula, modules=['numpy'])

    def __call__(self, *input_values):
        if isinstance(self.formula, bool) or (len(self.input_vars) == 0):
            return self.formula
        # print self.formula
        # print type(self.formula)
        if self.formula is not None:
            return self._lambdified_formula(*input_values)
        else:
            return self.boolean_outputs[sum(2**i * val for i, val in range(len(input_values)))]

    def __str__(self):
        return " " + str(self.formula)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, bool) or isinstance(other, sympy.boolalg.BooleanTrue) or \
           isinstance(other, sympy.boolalg.BooleanFalse) or (isinstance(other, int) and other in [0, 1]):
            if len(self.input_vars) == 0:
                return bool(self.formula) == bool(other)
            else:
                return False
        if isinstance(other, BooleanSymbolicFunc):
            return self.boolean_outputs == other.boolean_outputs
        try:
            for input_comb in itertools.product([False, True], repeat=len(self.input_vars)):
                if self(*input_comb) != other(*input_comb):
                    return False
        except ValueError:
            return False
        return True

    def __hash__(self):
        return hash(self.boolean_outputs)

    def __ne__(self, other):
        return not self == other

    def __nonzero__(self):
        if not isinstance(self.formula, bool) and not isinstance(self.formula, int) and \
           not isinstance(self.formula, (sympy.boolalg.BooleanTrue, sympy.boolalg.BooleanFalse)):
            raise ValueError("Cannot convert non constant BooleanSymbolicFunction to bool")
        if isinstance(self.formula, (sympy.boolalg.BooleanTrue, sympy.boolalg.BooleanFalse)):
            return self.formula == True
        return bool(self.formula)

    @staticmethod
    def sanitized_nand(*args):
        """
        Replaces sympy nand with Or of Nots (because Nand introduces problems with other replacements)
        :param args:
        :return:
        """
        return sympy.Or(*(sympy.Not(x) for x in args))

    def compose(self, input_funcs, simplify=True):
        """
        Composes symbolic boolean functions. Assumes input_funcs are ordered in the order of self.input_vars.
        After composition, returns the new function, with its inputs ordered by name.
        :param input_funcs:
        :param simplify:
        :return:
        """
        assert len(input_funcs) == len(self.input_vars)
        for f in input_funcs:
            if not isinstance(f, BooleanSymbolicFunc):
                raise NotImplementedError(
                    "Can't compose a symbolic boolean function with a function of type {}".format(f.type))

        print([f.formula for f in input_funcs])
        nand_free_formulas = [f.formula.replace(sympy.Nand,
                                                 BooleanSymbolicFunc.sanitized_nand) for f in input_funcs]

        replacement_dict = dict(zip(self.input_vars, nand_free_formulas))
        new_exp = self.formula.replace(sympy.Nand, BooleanSymbolicFunc.sanitized_nand).\
            subs(replacement_dict, simultaneous=True)
        if simplify:
            new_exp = sympy.simplify(new_exp)
        return BooleanSymbolicFunc(formula=new_exp)

    @staticmethod
    def from_sympy_func(sympy_func, variable_names):
        symbols = sympy.symbols(variable_names)
        expr = sympy_func(*symbols)
        return BooleanSymbolicFunc(formula=expr)


class SymmetricThresholdFunction(object):
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
        if other is None:
            return False
        elif other in [False, True, sympy.false, sympy.true]:
            other_func = lambda _: other
        else:
            other_func = other
        try:
            for input_comb in itertools.product([False, True], repeat=len(self.signs)):
                if self(*input_comb) != other_func(*input_comb):
                    return False
        except (ValueError, TypeError):
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


def perturb_line(f, line_indices, return_symbolic=False, n_inputs=None):
    """
    Given a logic function (possibly SymbolicBooleanFunction, but not necessarily), and an index of a truth
    table row to "perturb" (/flip), returns a function agreeing with the input function on all inputs but the line
    indices.
    If return_symbolic is true, creates a new SymbolicBooleanFunction. Otherwise just wraps the original one.
    :param f:
    :param line_indices:
    :param return_symbolic:
    :param n_inputs: if return_symbolic is true and f is not a symbolic boolean function,
    this specifies how many inputs f receives.
    :return:
    """

    if not return_symbolic:
        def perturbed_wrapper(*args):
            line_index = sum(2**i for i, b in enumerate(args) if b)
            return (1 - int(bool(f(*args)))) if line_index in line_indices else int(bool(f(*args)))
        return perturbed_wrapper
    else:
        if isinstance(f, BooleanSymbolicFunc):
            original_outputs = f.boolean_outputs
        else:
            original_outputs = [f(*args) for args in itertools.product([False, True], repeat=n_inputs)]
        boolean_outputs = list(original_outputs)
        for index in line_indices:
            boolean_outputs[index] = 1 - bool(boolean_outputs[index])  # to work with sympy's logic
        input_names = [x.name for x in f.input_vars] if isinstance(f, BooleanSymbolicFunc) else \
            ["var_{}".format(i) for i in range(n_inputs)]
        return BooleanSymbolicFunc(input_names=input_names,
                                               boolean_outputs=boolean_outputs)


def expression_without_variable(var_name, expression):
    """
    Removes any use of the variable name in the sympy expression. This is done by recursively removing the variable
    from every argument in the expression and its sub-expressions, and removing empty expressions that result.
    If the entire expression is empty, returns None (note this also converts sympy.false and sympy.true to None)
    :param expression:
    :return:
    """
    if expression.is_symbol:
        return expression if expression.name != var_name else None
    new_args = [expression_without_variable(var_name, arg) for arg in expression.args]
    new_args = [arg for arg in new_args if arg is not None]
    if len(new_args) == 0:
        return None
    return expression.func(*new_args)