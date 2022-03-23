#!/usr/bin/env python3

from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from ladder_network_filter import LadderNetworkFilter

def polynomial_array_strip_high_order_zeros(x,tol=1e-10):
    while len(x) > 0:
        if abs(x[-1]) > tol:
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
    from scipy import signal
    from filter_design import semi_gaussian_complex_all_pole_filter, semi_gaussian_complex_pole_locations, plot_filters_behavior

    real_pole_filter2 = cauerI_synthesis([1],[1,2,1])
    real_pole_filter3 = cauerI_synthesis([1],[8,12,6,1])
    real_pole_filter4 = cauerI_synthesis([1],[81,108,54,12,1])
    LadderNetworkFilter.make_plots_many_filters([real_pole_filter2,real_pole_filter3,real_pole_filter4],[r"$\frac{1}{(s+1)^2}$",r"$\frac{1}{(s+2)^3}$",r"$\frac{1}{(s+3)^4}$"],"Real_Pole_Filters.pdf","Cauer I LC Filters",1e-3,1e3,1e-3,0,5)

    bessel_filters = []
    for i in range(2,10):
        n, d = signal.bessel(i,1,btype="lowpass",analog=True,output="ba")
        #p, z, k = signal.tf2zpk(n,d)
        #p /= i-1
        #z /= i-1
        #n, d = signal.zpk2tf(p,z,k)
        ladder_filter = cauerI_synthesis(n,d)
        print(ladder_filter)
        bessel_filters.append(ladder_filter)
    bessel_titles = ["Bessel {}O".format(i) for i in range(2,10)]
    LadderNetworkFilter.make_plots_many_filters(bessel_filters,bessel_titles,"Bessel_Filters.pdf","Cauer I LC Bessel Filters",1e-3,1e3,1e-3,0,5)
    
    semi_gaus_filters = []
    semi_gaus_titles = ["Semi-Gaus {}O".format(i) for i in range(3,6)]
    for i in range(3,6):
        poles = semi_gaussian_complex_pole_locations(i)/(i-2)/4
        zpg = signal.ZerosPolesGain([],poles,[1])
        tf = zpg.to_tf()
        ladder_filter = cauerI_synthesis(tf.num,tf.den)
        semi_gaus_filters.append(ladder_filter)
    LadderNetworkFilter.make_plots_many_filters(semi_gaus_filters,semi_gaus_titles,"Synth_Semi_Gaus.pdf","Cauer I LC Filters",1e-3,1e3,1e-3,0,7)

    LadderNetworkFilter.make_plots_many_filters([bessel_filters[1],semi_gaus_filters[0],real_pole_filter3],[bessel_titles[1],semi_gaus_titles[0],r"$\frac{1}{(s+2)^3}$"],"Synth_Comparison.pdf","3rd Order LC Filter Comparison",1e-3,1e3,1e-3,0,7)
