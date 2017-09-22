import math

def list_repr(elements):
    if len(elements) == 0:
        return "[]"
    else:
        return "[" + str(reduce(lambda x, y: str(x) + ", " + str(y), elements)) + "]"


def slice_int(x, k):
    """
    Slices an integer to k-bit parts, and returns the respective parts (in range 0-2**k - 1), MSB to LSB.
    :param x: A positive integer
    :param k: An integer representing the size of the slices, in bits
    :return:
    """
    return [(x >> i*k) % 2**((i + 1)*k) for i in range(int(math.ceil(math.log(x, 2)/k)))][::-1]

