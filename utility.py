import fractions
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