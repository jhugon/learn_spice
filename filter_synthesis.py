#!/usr/bin/env python3

from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

def polynomial_array_strip_high_order_zeros(x):
    while len(x) > 0:
        if x[-1] != 0.:
            break
        x = x[:-1]
    return x

def polynomial_divide(n,d):
    """
    Performs polynomial division

    n: an array of numerator coefficients, where n[i] corresponds to the
        coefficient of the s^i term

    d: an array of denominator coefficients, where d[i] corresponds to the
        coefficient of the s^i term

    returns (q,rn) where q is the quotient array of coefficients, rn is the
        remainder numerator (it can be divided by d to get the full remainder)
    """

    n = np.array(n,dtype="float64")
    d = np.array(d,dtype="float64")
    # strip higher-order terms that are zero
    n = polynomial_array_strip_high_order_zeros(n)
    d = polynomial_array_strip_high_order_zeros(d)
    assert(len(n)>=len(d))

    qscaler = n[-1]/d[-1]
    qorder = len(n)-len(d)
    q = np.zeros(qorder+1)
    q[-1] = qscaler

    # now multiply q * d and subtract from n to find the remainder
    mult = np.zeros(len(d)+qorder)
    mult[qorder:] = qscaler * d
    rn = n-mult
    # strip higher-order terms that are zero
    q = polynomial_array_strip_high_order_zeros(q)
    rn = polynomial_array_strip_high_order_zeros(rn)
    return q, rn

def polynomial_continued_fraction_decomp(n,d):
    """
    Performs continued fraction decomposition of the fraction represented by
        numerator n and denominator d

    n: an array of numerator coefficients, where n[i] corresponds to the
        coefficient of the s^i term

    d: an array of denominator coefficients, where d[i] corresponds to the
        coefficient of the s^i term

    returns a list containing the coefficients of the continued fraction. Each
        entry in the list represents a polynomial, where the ith entry corresponds to
        the coefficient of the s^i term.
    """

    q, rn = polynomial_divide(n,d)
    result = [q]
    if len(rn) == 0: # no remainder
        return result
    else:
        return result + polynomial_continued_fraction_decomp(d,rn)

if __name__ == "__main__":
    n = [0,24,0,30,0,9]
    d = [8,0,36,0,18]
    print(polynomial_divide(n,d))
    print(polynomial_continued_fraction_decomp(n,d))
