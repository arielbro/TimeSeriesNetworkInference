import math


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

