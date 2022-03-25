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

def polynomial_multiply(x1,x2):
    """
    Performs polynomial multiplication

    x1, x2: arrays of coefficients, where x[i] corresponds to the
        coefficient of the s^i term

    result is a similar array of coefficients to x1 and x2
    """

    x1 = np.array(x1,dtype="float64")
    x2 = np.array(x2,dtype="float64")
    # strip higher-order terms that are zero
    x1 = polynomial_array_strip_high_order_zeros(x1)
    x2 = polynomial_array_strip_high_order_zeros(x2)

    if len(x1) == 0:
        return x1
    if len(x2) == 0:
        return x2

    len1 = len(x1)
    len2 = len(x2)
    order1 = len1-1
    order2 = len2-1
    order_res = order1+order2
    len_res = order_res+1
    tmp = np.zeros((len1,len_res),dtype="float64")
    for shift1 in range(len1):
        tmp[shift1,shift1:shift1+len2] = x1[shift1]*x2
    result = tmp.sum(axis=0)
    return result


def polynomial_divide(n,d):
    """
    Performs polynomial division.

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
    order_n = len(n)-1
    order_d = len(d)-1
    quotient_order = order_n-order_d
    quotient_len = quotient_order+1
    if quotient_len <= 0:
        return np.zeros(1,dtype="float64"), n
    quotient = np.zeros(quotient_len,dtype="float64")

    remainder_n = np.array(n)
    for iTerm in reversed(range(quotient_len)):
        q_term, remainder_n = polynomial_divide_just_lead_term(remainder_n,d)
        quotient[:iTerm+1] = q_term
    return quotient, remainder_n



def polynomial_divide_just_lead_term(n,d):
    """
    Performs part of polynomial division. The quotient is only the leading term
        of the quotient, and the remainder is the remainder of that.

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
    if len(n) < len(d):
        return np.zeros(1,dtype="float64"), n

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

def polynomial_sqrt(x):
    """
    Performs polynomial square root.

    x: an array of numerator coefficients, where n[i] corresponds to the
        coefficient of the s^i term
    """

    x = np.array(x,dtype="complex128")
    # strip higher-order terms that are zero
    x = polynomial_array_strip_high_order_zeros(x)
    len_x = len(x)
    order_x = len_x -1
    if order_x % 2 == 1:
        raise ValueError(f"Polynomial must have even order: {x}")
    result_order = int(order_x/2)
    assert(result_order*2 == order_x)
    result_len = result_order+1
    result = np.zeros(result_len,dtype="complex128")
    result[-1] = np.sqrt(x[-1])
        

def cauerI_synthesis_equal_inout_impedance(n,d,R=1.,reverse_polys=False):
    """
    Use Cauer I synthesis to make a LC low-pass filter from the Vout/Vin
        transfer function given by the numerator, n, and denominator, d, polynomials.
        Constructs the ladder shunt (capacitor) first.

    Williams, A. Analog Filter and Circuit Design Handbook. McGraw Hill (2013).
        ISBN-13: 978-0071816717
        Section 1.2.1

    Assumes the numerator is just [1]

    n: an array of numerator coefficients, where n[i] corresponds to the
        coefficient of the s^i term. If reverse_polys is True, then the
        coefficients are in the opposite order.

    d: an array of denominator coefficients, where d[i] corresponds to the
        coefficient of the s^i term. If reverse_polys is True, then the
        coefficients are in the opposite order.

    R: the assumed source and load impedance

    returns (C,L), where C is an array of capacitances in F and L is an array
        of inductances in H.
    """

    n = np.array(n,dtype="float64")
    d = np.array(d,dtype="float64")
    if reverse_polys:
        n = np.flip(n)
        d = np.flip(d)
    # strip higher-order terms that are zero
    n = polynomial_array_strip_high_order_zeros(n)
    d = polynomial_array_strip_high_order_zeros(d)
    assert(len(n)<=len(d))

    # Only works for just 1 in the numerator for now
    assert(len(n) == 1)
    assert(n[0] == 1.)


    d_squared = polynomial_multiply(d,d)
    K_squared_n = -np.array(d_squared)
    K_squared_n[0] += 1.
    K_squared_d = d_squared
    breakpoint()

    driving_point_Z_n = np.array(d)
    driving_point_Z_d = np.array(d)

    driving_point_Z_n[-1] -= 1.
    driving_point_Z_d[-1] += 1.

    cfd = polynomial_continued_fraction_decomp(driving_point_Z_n,driving_point_Z_d)
    cfd_y = polynomial_continued_fraction_decomp(driving_point_Z_d,driving_point_Z_n)

    if len(cfd[0]) == 1 and abs(cfd[0][0]) < 1e-6:
        cfd.pop(0)
    else:
        raise Exception(f"First element of continued fraction should be [0], not: {cfd[0]}")
    if len(cfd[-1]) == 1 and abs(cfd[-1][0]-1.) < 1e-6:
        cfd.pop(-1)
    else:
        raise Exception(f"Last element of continued fraction should be [1], not: {cfd[-1]}")
    for x in cfd:
        assert(len(x) == 2)

    cfd_s_coefs = np.array([x[1] for x in cfd],dtype="float64")

    C = cfd_s_coefs[::2]
    L = cfd_s_coefs[1::2]
    C /= R
    L *= R

    lnf = LadderNetworkFilter(C,L,Rin=R,Rout=R,shunt_first=True)
    return lnf


def cauerI_synthesis_inf_in_impedance(n,d,R=1.,reverse_polys=False):
    """
    Use Cauer I synthesis to make a LC low-pass filter from the trans-impedance
        (Z_{21}) transfer function given by the numerator, n, and denominator, d,
        polynomials. Since the trans-impedance is converted, the filter will be
        shunt-capacitor first.

    **These filters are designed for infinite input impedance**

    Bakshi, U. A., Bakshi, A. V. Network Analysis & Synthesis: Laplace Transform, 
        Two Port Networks, Network Synthesis. Technical Publications (2020).
        ISBN-13: 978-9333223515

    Assumes the numerator is just [1]

    n: an array of numerator coefficients, where n[i] corresponds to the
        coefficient of the s^i term. If reverse_polys is True, then the
        coefficients are in the opposite order.

    d: an array of denominator coefficients, where d[i] corresponds to the
        coefficient of the s^i term. If reverse_polys is True, then the
        coefficients are in the opposite order.

    R: the assumed source and load impedance

    returns (C,L), where C is an array of capacitances in F and L is an array
        of inductances in H.
    """

    n = np.array(n,dtype="float64")
    d = np.array(d,dtype="float64")
    if reverse_polys:
        n = np.flip(n)
        d = np.flip(d)
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
    cfd_s_coefs = np.array([x[1] for x in cfd],dtype="float64")
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
    import sys

    x = polynomial_multiply([1,2,3],[3,4])
    y = polynomial_divide(x,[3,4])
    z = polynomial_divide(x,[1,2,3])
    cf = polynomial_continued_fraction_decomp([0,24,0,30,0,9],[8,0,36,0,18])
    f = cauerI_synthesis_inf_in_impedance([1],[2,3,3,1])
    f = cauerI_synthesis_equal_inout_impedance([1],[2,3,3,1])
    sys.exit()

    n = 3
    print(signal.butter(n,1,btype="lowpass",analog=True,output="ba"))
    print(signal.bessel(n,1,btype="lowpass",analog=True,output="ba"))
    print(signal.cheby1(n,0.1,1,btype="lowpass",analog=True,output="ba"))
    print(signal.cheby2(n,10,1,btype="lowpass",analog=True,output="ba"))


    # 2nd order 0.1dB Chebyshev
    zpg = signal.ZerosPolesGain([],[-0.6743+0.7075j,-0.6743-0.7075j],[1])
    print(zpg)
    tf = zpg.to_tf()
    print(tf)
    ladder_filter = cauerI_synthesis_equal_inout_impedance(tf.num,tf.den,reverse_polys=True)
    print(ladder_filter)

    def get_butterworth(n):
        return signal.butter(n,1,btype="lowpass",analog=True,output="ba")
    def get_bessel(n):
        return signal.bessel(n,1,btype="lowpass",analog=True,output="ba")

    sys.exit(1)

    for i in range(2,8):
        ladder_filter = cauerI_synthesis_equal_inout_impedance(*get_bessel(i),reverse_polys=True)
        print(ladder_filter)

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
