#!/usr/bin/env python3

from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from ladder_network_filter import LadderNetworkFilter

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

def cauerI_synthesis(n,d,R=1.):
    """
    Use Cauer I synthesis to make a LC low-pass filter from the trans-impedance
        (Z_{21}) transfer function given by the numerator, n, and denominator, d,
        polynomials. Since the trans-impedance is converted, the filter will be
        shunt-capacitor first.

    Assumes the numerator is just [1]

    n: an array of numerator coefficients, where n[i] corresponds to the
        coefficient of the s^i term

    d: an array of denominator coefficients, where d[i] corresponds to the
        coefficient of the s^i term

    R: the assumed source and load impedance

    returns (C,L), where C is an array of capacitances in F and L is an array
        of inductances in H.
    """

    n = np.array(n,dtype="float64")
    d = np.array(d,dtype="float64")
    # strip higher-order terms that are zero
    n = polynomial_array_strip_high_order_zeros(n)
    d = polynomial_array_strip_high_order_zeros(d)
    assert(len(n)<=len(d))

    # Only works for just 1 in the numerator for now
    assert(len(n) == 1)
    assert(n[0] == 1.)

    range_mod_2 = np.arange(len(d)) % 2
    d_evens = polynomial_array_strip_high_order_zeros(d * (range_mod_2 == 0))
    d_odds = polynomial_array_strip_high_order_zeros(d * (range_mod_2 == 1))
    zfirst = len(d_evens) >= len(d_odds) # first component is impedance, otherwise is admittance (y)

    cfd = []
    if zfirst:
        cfd = polynomial_continued_fraction_decomp(d_evens,d_odds)
    else:
        cfd = polynomial_continued_fraction_decomp(d_odds,d_evens)

    for x in cfd:
        assert(len(x) == 2)
    cfd_s_coefs = np.array([x[1] for x in cfd])
    term_is_z = (np.arange(len(cfd_s_coefs)) % 2)
    if zfirst:
        term_is_z = np.logical_not(term_is_z)

    # cfd_s_coefs and term_is_z both go from output to input port, so reverse
    terms_flipped = np.flip(cfd_s_coefs)
    term_is_z_flipped = np.flip(term_is_z)
    # Closest element to input port is capacitor, so first (and even) elements are caps
    C = terms_flipped[::2]
    L = terms_flipped[1::2]
    if term_is_z_flipped[0]:
        C = 1./C
    if not term_is_z_flipped[1]:
        L = 1./L
    C /= R
    L *= R
    lnf = LadderNetworkFilter(C,L,Rin=R,Rout=R,shunt_first=True)
    return lnf

if __name__ == "__main__":
    real_pole_filter2 = cauerI_synthesis([1],[1,2,1])
    print("1/(s+1)^2 filter: ",real_pole_filter2)
    print("Poles: {}, Zeros: {}".format(*real_pole_filter2.get_poles_zeros()))
    #real_pole_filter2.make_plots("Real_Pole_Filter2.pdf",r"$\frac{1}{(s+1)^2}$ Cauer I Synthesis",1e-3,1e3,1e-3,0,5)
    real_pole_filter3 = cauerI_synthesis([1],[8,12,6,1])
    print("1/(s+2)^3 filter: ",real_pole_filter3)
    print("Poles: {}, Zeros: {}".format(*real_pole_filter3.get_poles_zeros()))
    #real_pole_filter3.make_plots("Real_Pole_Filter3.pdf",r"$\frac{1}{(s+2)^3}$ Cauer I Synthesis",1e-3,1e3,1e-3,0,5)
    real_pole_filter4 = cauerI_synthesis([1],[81,108,54,12,1])
    print("1/(s+3)^4 filter: ",real_pole_filter4)
    print("Poles: {}, Zeros: {}".format(*real_pole_filter3.get_poles_zeros()))
    #real_pole_filter3.make_plots("Real_Pole_Filter4.pdf",r"$\frac{1}{(s+3)^4}$ Cauer I Synthesis",1e-3,1e3,1e-3,0,5)

    LadderNetworkFilter.make_plots_many_filters([real_pole_filter2,real_pole_filter3,real_pole_filter4],[r"$\frac{1}{(s+1)^2}$",r"$\frac{1}{(s+2)^3}$",r"$\frac{1}{(s+3)^4}$"],"Real_Pole_Filters.pdf","Cauer I LC Filters",1e-3,1e3,1e-3,0,5)
