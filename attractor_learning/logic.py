import sympy
import itertools
import math
from .utility import list_repr


class BooleanSymbolicFunc(object):
    # Printing lists the truth table only up to this many inputs; beyond it, just the count of True rows.
    _MAX_PRINTED_INPUTS = 6

    def __init__(self, input_names=None, boolean_outputs=None):
        # make all fields immutable, so the function can be shallow copied safely.
        # The truth table is the function: evaluation is a table lookup, first input as the most significant
        # bit (see __call__).
        boolean_outputs = tuple(boolean_outputs)
        if input_names is None:
            n_inputs = int(math.log(len(boolean_outputs), 2))
            input_names = ['input_{}'.format(i) for i in range(n_inputs)]

        if len(input_names) != math.frexp(len(boolean_outputs))[1] - 1:
            raise ValueError("non-matching length for variable names list and boolean outputs list")
        # assumes boolean_outputs is a power of 2
        self.input_vars = tuple(sympy.symbols(name) for name in input_names)
        self._boolean_outputs = boolean_outputs

    @property
    def boolean_outputs(self):
        return self._boolean_outputs

    def __call__(self, *input_values):
        # Straight truth-table lookup, first input as the most significant bit - the order boolean_outputs
        # is built in (itertools.product over the inputs) everywhere it is read. With no inputs the table is
        # the one constant row.
        row = 0
        for value in input_values:
            row = (row << 1) | (1 if value else 0)
        return self._boolean_outputs[row]

    def __str__(self):
        names = ", ".join(var.name for var in self.input_vars)
        if len(self.input_vars) <= BooleanSymbolicFunc._MAX_PRINTED_INPUTS:
            table = "".join("1" if out else "0" for out in self._boolean_outputs)
            return " f({}) = {}".format(names, table)
        return " f({}) = {} of {} rows True".format(names, sum(1 for out in self._boolean_outputs if out),
                                                    len(self._boolean_outputs))

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, bool) or isinstance(other, sympy.boolalg.BooleanTrue) or \
           isinstance(other, sympy.boolalg.BooleanFalse) or (isinstance(other, int) and other in [0, 1]):
            if len(self.input_vars) == 0:
                return bool(self._boolean_outputs[0]) == bool(other)
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

    def to_dict(self):
        """Serialization via the (input_names, boolean_outputs) form, which round-trips exactly through the
        constructor and preserves input order. These functions are always low in-degree in this codebase, so
        the truth table is cheap (unlike threshold functions, which can have a high in-degree)."""
        return {"type": "boolean_symbolic",
                "input_names": [var.name for var in self.input_vars],
                "boolean_outputs": [bool(out) for out in self.boolean_outputs]}

    @staticmethod
    def from_dict(d):
        return BooleanSymbolicFunc(input_names=d["input_names"], boolean_outputs=d["boolean_outputs"])


class SparseBooleanFunc(object):
    """A general Boolean function stored as the truth-table rows that disagree with its majority output.

    Fills the same role as BooleanSymbolicFunc - an arbitrary Boolean function of named inputs, called
    positionally in predecessor order - but does not hold a 2**in-degree truth table.
    A strongly biased function, which is what criticality demands at a high in-degree, disagrees with its
    majority value on few rows, so storing those row indices costs O(minority rows) rather than
    O(2**in-degree), and neither construction nor evaluation ever walks the full table.

    `minority_rows` index the truth table with the FIRST input as the most significant bit - the order
    BooleanSymbolicFunc.boolean_outputs is built in - so the two representations number rows identically
    and a function can be converted between them without reordering.
    """

    def __init__(self, input_names, minority_rows, default_output):
        self.input_names = tuple(input_names)
        self.default_output = bool(default_output)
        self.minority_rows = frozenset(int(row) for row in minority_rows)
        n_rows = 2 ** len(self.input_names)
        if self.minority_rows and (min(self.minority_rows) < 0 or max(self.minority_rows) >= n_rows):
            raise ValueError("minority row index out of range for {} inputs".format(len(self.input_names)))
        # The stored form is canonical, and that is what lets __eq__ and __hash__ answer without
        # materializing 2**in-degree rows: once the minority is required to be the smaller side of the
        # table, with True taking the tie at exactly half, a function has exactly one description here.
        if (2 * len(self.minority_rows) > n_rows) or \
                (2 * len(self.minority_rows) == n_rows and not self.default_output):
            raise ValueError("minority_rows must be the smaller side of the table (got {} of {} rows with "
                             "default {}); a function that balanced belongs in a BooleanSymbolicFunc"
                             .format(len(self.minority_rows), n_rows, self.default_output))

    def __call__(self, *input_values):
        if len(input_values) != len(self.input_names):
            raise ValueError("expected {} inputs, got {}".format(len(self.input_names), len(input_values)))
        row = 0
        for value in input_values:
            row = (row << 1) | (1 if value else 0)
        return (not self.default_output) if row in self.minority_rows else self.default_output

    @property
    def boolean_outputs(self):
        """The full truth table in BooleanSymbolicFunc's row order. 2**in-degree entries, which is exactly
        what this class exists to avoid - it is here for the small-in-degree comparisons that want it."""
        return tuple((not self.default_output) if row in self.minority_rows else self.default_output
                     for row in range(2 ** len(self.input_names)))

    def _canonical_key(self):
        """A hashable key identical iff two functions are equal - see the canonical form in __init__."""
        return (len(self.input_names), self.default_output, self.minority_rows)

    def to_dict(self):
        """Serialization by minority rows; sorted so the file is byte-identical across runs (the rows are
        held in a set, whose iteration order is not a promise)."""
        return {"type": "sparse_boolean",
                "input_names": list(self.input_names),
                "default_output": self.default_output,
                "minority_rows": sorted(self.minority_rows)}

    @staticmethod
    def from_dict(d):
        return SparseBooleanFunc(input_names=d["input_names"], minority_rows=d["minority_rows"],
                                 default_output=d["default_output"])

    def __str__(self):
        return "default={}, {} minority row(s) of {}".format(
            self.default_output, len(self.minority_rows), 2 ** len(self.input_names))

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if other is None:
            return False
        # compact comparisons (no 2**in-degree truth table) for the cases that occur in bulk
        if isinstance(other, SparseBooleanFunc):
            return self._canonical_key() == other._canonical_key()
        if other in [False, True, sympy.false, sympy.true]:
            return len(self.minority_rows) == 0 and self.default_output == bool(other)
        # general fallback for arbitrary callables, at the cost this class is built to avoid
        try:
            for input_comb in itertools.product([False, True], repeat=len(self.input_names)):
                if self(*input_comb) != other(*input_comb):
                    return False
        except (ValueError, TypeError):
            return False
        return True

    def __hash__(self):
        return hash(self._canonical_key())

    def __ne__(self, other):
        return not self == other


class SymmetricThresholdFunction(object):
    # TODO: implement in ILP model finding (threshold is not boolean, not supported there ATM)
    def __init__(self, signs, threshold):
        # translate signs to bool values, if not already. Signs must be nonzero (+1/-1); a zero sign means
        # "input unused" and is only meaningful inside the inference ILP, not in a realized function.
        self.signs = []
        for sign in signs:
            if isinstance(sign, bool):
                self.signs.append(sign)
            elif isinstance(sign, int):
                assert sign in [1, -1], "signs must be nonzero (+1/-1), got {}".format(sign)
                self.signs.append(True if sign == 1 else False)
            else:
                raise ValueError("illegal type for signs:{}".format(type(sign)))
        # threshold in [0, degree+1]: 0 is constant-True, degree+1 is constant-False, [1, degree] is non-constant.
        assert 0 <= threshold <= len(self.signs) + 1, \
            "threshold {} out of range [0, {}]".format(threshold, len(self.signs) + 1)
        self.threshold = threshold

    def _canonical_key(self):
        """A hashable key identical iff two functions are equal. With nonzero signs, a non-constant function
        (1 <= threshold <= degree) has a unique (signs, threshold); the constant cases are independent of the
        signs, so collapse threshold <= 0 (always True) and threshold > degree (always False)."""
        degree = len(self.signs)
        if self.threshold <= 0:
            return ("const", True)
        if self.threshold > degree:
            return ("const", False)
        return ("threshold", tuple(self.signs), self.threshold)

    def to_dict(self):
        """Portable, truth-table-free serialization."""
        return {"type": "symmetric_threshold",
                "signs": [1 if sign else -1 for sign in self.signs],
                "threshold": self.threshold}

    @staticmethod
    def from_dict(d):
        return SymmetricThresholdFunction(signs=d["signs"], threshold=d["threshold"])

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
        # compact comparisons (no 2**degree truth table) for the common cases
        if isinstance(other, SymmetricThresholdFunction):
            return self._canonical_key() == other._canonical_key()
        if other in [False, True, sympy.false, sympy.true]:
            key = self._canonical_key()
            return key[0] == "const" and key[1] == bool(other)
        # general fallback for arbitrary callables (rare; only reached for non-threshold function types)
        try:
            for input_comb in itertools.product([False, True], repeat=len(self.signs)):
                if self(*input_comb) != other(*input_comb):
                    return False
        except (ValueError, TypeError):
            return False
        return True

    def __hash__(self):
        return hash(self._canonical_key())

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
