import numpy
import fractions
import math
import sympy
from functools import reduce


def rotate(l, n):
    return l[n:] + l[:n]


class Attractor:
    # TODO: Use a general class for attractors (everywhere)
    def __init__(self, states):
        largest_ind = max(list(range(len(states))), key=lambda t: order_key_func(states[t]))
        self.states = tuple(tuple(states[(t + largest_ind + 1) % len(states)]) for t in range(len(states)))

    def __eq__(self, other):
        """
        Compares attractors invariant to rotation (by rotating both to have largest state last)
        Note that the container types used for attractor states are important (e.g. [0,0,0] != (0,0,0))
        :param self:
        :param other:
        :return:
        """
        if not isinstance(other, Attractor):
            raise NotImplementedError("Can't compare an attractor to anything else")
        if len(self.states) != len(other.states):
            return False
        return self.states == other.states

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return tuple([tuple(s) for s in self.states]).__hash__()

    def __str__(self):
        return "Attractor, states: {}".format(self.states)

    def __repr__(self):
        return str(self)

def order_key_func(node_states): return sum(node * 2 ** i for (i, node) in enumerate(node_states))


def list_repr(elements):
    if len(elements) == 0:
        return "[]"
    else:
        return "[" + str(reduce(lambda x, y: str(x) + ", " + str(y), elements)) + "]"


def slice_int(x, k, d):
    """
    Slices an integer in the domain [0,2**d] to k-bit parts, and returns the
    respective parts (in range 0-2**k - 1), MSB to LSB.
    :param x: A positive integer
    :param k: An integer representing the size of the slices, in bits
    :param d: The range of x. If x < 2**d - 1, it is treated as left padded with zeros.
    :return:
    """
    return [(x >> i*k) % 2**k for i in range(int(math.ceil(float(d)/k)))][::-1]


def divisors(x):
    """
    divisor list for n, copy-pasted (and corrected for 1) from stack-overflow
    :param x:
    :return:
    """
    if x == 1:
        return [x]
    div_list = []
    y = 1
    while y <= int(math.sqrt(x)):
        if x % y == 0:
            div_list.append(y)
            if x / y != y:
                div_list.append(int(x / y))
        y += 1
    return div_list


def phi(n):
    """
    euler's totient function, copy-pasted from stack-overflow
    :param n:
    :return:
    """
    amount = 0
    for k in range(1, n + 1):  # n+1 important for phi(1)=1
        if fractions.gcd(n, k) == 1:
            amount += 1
    return amount


def binary_necklaces(n):
    """
    returns the number of binary necklaces over n bits, i.e. number of binary strings up to rotation.
    derivation here - www.quora.com/How-many-unique-binary-matrices-are-there-up-to-rotations-translations-and-flips
    also here - https://en.wikipedia.org/wiki/Necklace_(combinatorics)#Number_of_necklaces
    :param n:
    :return:
    """
    s = 0
    for divisor in divisors(n):
        s += phi(n / divisor) * (2**divisor)
        # print "d={}, phi(n/d)={}, s={}".format(divisor, phi(n/divisor), s)
    return s / n


def attractor_sets_equality(first_attractors, second_attractors):
    """
    Compares two containers of attractors, invariant to rotations.
    For efficiency (hopefully), creates a set-like structure for comparisons.
    :param self.statess:
    :param other.statess:
    :return:
    """
    if len(first_attractors) != len(second_attractors):
        return False
    first_attractors_set = set(Attractor(att) for att in first_attractors)
    second_attractors_set = set(Attractor(att) for att in second_attractors)
    return first_attractors_set == second_attractors_set


def is_attractor_in_attractor_list(attractor, attractor_list):
    first_attractor = Attractor(attractor)
    second_attractors_set = set(Attractor(att) for att in attractor_list)
    return first_attractor in second_attractors_set


def is_same_attractor(a1, a2):
    """
    :param a1: an attractor, represented as an ordered iterable of network states
    :param a2: ""
    :return: True iff the attractors have same states in same order, up to a shift.
    """
    if len(a1) != len(a2):
        return False
    a1 = tuple(tuple(1 if v_state else 0 for v_state in s) for s in a1)
    a2 = tuple(tuple(1 if v_state else 0 for v_state in s) for s in a2)
    for shift in range(len(a1)):
        if a1 == rotate(a2, shift):
            return True
    return False


def is_same_state(s1, s2):
    """
    :param s1: A network state, represented as an ordered iterable of values interpretable as boolean.
    :param s2: ""
    :return: True if s1 and s2 represent the same state
    """
    if len(s1) != len(s2):
        # return False
        raise ValueError("Can't compare states from models of different size.")
    for v_state in s1:
        assert v_state in [0, 1, False, True, sympy.false, sympy.true]
    for v_state in s2:
        assert v_state in [0, 1, False, True, sympy.false, sympy.true]
    s1_standard = tuple(1 if v_state else 0 for v_state in s1)
    s2_standard = tuple(1 if v_state else 0 for v_state in s2)
    return s1_standard == s2_standard


def is_attractor_valid(attractor, G):
    """
    Checks whether an attractor is valid in the model G, regardless of representation.
    :param attractor:
    :param G:
    :return:
    """
    # TODO: write tests!
    for state, next_state in zip(attractor, rotate(attractor, 1)):
        if not is_same_state(G.next_state(state), next_state):
            return False
    return True
