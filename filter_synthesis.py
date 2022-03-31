#!/usr/bin/env python3

from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import Polynomial
from ladder_network_filter import LadderNetworkFilter

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

    if not isinstance(n,Polynomial):
        n = Polynomial(n)
    if not isinstance(d,Polynomial):
        d = Polynomial(d)
    q, rn = divmod(n,d)
    q = q.trim()
    rn = rn.trim()
    result = [q]
    if len(rn) == 0 or (len(rn) == 1 and abs(rn.coef[0]) < 1e-12): # no remainder
        return result
    else:
        return result + polynomial_continued_fraction_decomp(d,rn)

def polynomial_conj(p):
    """
    Assumes polynomial variable, s, is complex

    Input polynomial is an instance of np.polynomial.Polynomial

    returns an instance of np.polynomial.Polynomial
    """
    assert(isinstance(p,Polynomial))
    roots = p.roots()
    if p.degree() == 0:
        return p
    result = Polynomial.fromroots(-roots)
    return result

def polynomial_complex_sqrt(p):
    """
    Performs something like a polynomial square root.
    Finds f such that f times conj(f) = p

    Keep in mind -result is also a solution

    p: an instance of np.polynomial.Polynomial

    returns an instance of np.polynomial.Polynomial
    """

    tol = 1e-12

    p = p.trim()
    p = Polynomial([(coef if abs(coef) > tol else 0.) for coef in p.coef])
    roots = p.roots()
    roots = np.sort(roots)
    if len(roots) != p.degree():
        raise Exception(f"Found a different number of roots than expected. Roots: {roots} for polynomial: {p}")
    if np.all(roots < tol): # just zeros
        new_degree = p.degree() // 2
        new_coef = np.sqrt(abs(p.coef[-1]))
        new_coefs = [0.]*new_degree + [new_coef]
        result = Polynomial(new_coefs)
    else:
        roots_real_gt0 = roots[roots.real > 0]
        roots_real_lt0 = roots[roots.real < 0]
        assert(len(roots_real_gt0) == len(roots_real_lt0)) # assume real >0 roots are images of (good) neg roots

        result_roots = np.array(roots_real_lt0)
        for i in range(len(result_roots)):
            if abs(result_roots[i].imag) < tol:
                result_roots[i] = -abs(result_roots[i])
        result_poly = Polynomial.fromroots(result_roots)
        result_coef = result_poly.coef
        for i in range(len(result_coef)):
            if abs(result_coef[i].imag) < tol:
                result_coef[i] = -abs(result_coef[i])
        if abs(result_coef.imag).sum() < tol:
            result_coef = result_coef.real
        result = Polynomial(result_coef)

    result_conj = polynomial_conj(result)
    check_pos = result*result_conj
    check_neg = -result*result_conj
    error_pos = abs(p.coef-check_pos.coef)
    error_neg = abs(p.coef-check_neg.coef)
    if error_pos.sum() > tol and error_neg.sum() > tol:
        raise Exception(f"+/- result times it's conjugate doesn't match input: {p} result*conj(result): {check_pos} result: {result}")
    return result

def cauerI_synthesis_equal_inout_impedance(n,d,R=1.,reverse_polys=False):
    """
    Use Cauer I synthesis to make a LC low-pass filter from the Vout/Vin
        transfer function given by the numerator, n, and denominator, d, polynomials.
        Constructs the ladder shunt (capacitor) first.

    Youla, Dante C. Theory and Synthesis of Linear Passive Time-Invariant
        Networks. Cambridge University Press, 2015. 
        DOI: https://doi.org/10.1017/CBO9781316403105

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

    tol = 1e-12

    def only_linear_term_of_cfd(cfd):
        """
        returns array of linear term coefs and the final term's constant coef
        """
        result = np.zeros(len(cfd))
        result_const = None
        for iPoly, poly in enumerate(cfd):
            if len(poly) == 1:
                if abs(poly.coef[0]) < tol:
                    result[iPoly] = 0.
                    continue
                else:
                    raise Exception(f"Term {poly} in cfd: {cfd} is order 0 and the constant isn't zero.")
            if len(poly) != 2:
                raise Exception(f"Term {poly} in cfd: {cfd} has the wrong order, should be just 1.")
            if iPoly == len(cfd)-1:
                result_const = poly.coef[0] if abs(poly.coef[0]) > tol else 0.
            elif abs(poly.coef[0]) > tol:
                raise Exception(f"Term {poly} in cfd: {cfd} has a large constant coefficient")
            result[iPoly] = poly.coef[1]
        return result, result_const

    n = np.array(n,dtype="float64")
    d = np.array(d,dtype="float64")
    if reverse_polys:
        n = np.flip(n)
        d = np.flip(d)
    # strip higher-order terms that are zero
    n = Polynomial(n).trim()
    d = Polynomial(d).trim()
    assert(len(n)<=len(d))

    # Only works for just 1 in the numerator for now
    assert(len(n) == 1)
    assert(n.coef[0] == 1.)

    # power transfer function
    Gamma_n = n*polynomial_conj(n)
    Gamma_d = d*polynomial_conj(d)

    S_squared_n = -Gamma_d - Gamma_n
    S_squared_d = -Gamma_d

    # Now need to find S, where S times conj(S) = S_squared
    # Make sure to use +/- S
    S_n = polynomial_complex_sqrt(S_squared_n)
    S_d = polynomial_complex_sqrt(S_squared_d)

    # Driving point Z = (1+S)/(1-S)
    # That may be specifically for 1 ohm in/out impedance
    driving_point_Z_n_plus = S_d + S_n
    driving_point_Z_d_plus = S_d - S_n
    ## Use both + and - solution of polynomial_complex_sqrt
    driving_point_Z_n_minus = S_d - S_n
    driving_point_Z_d_minus = S_d + S_n

    cfd_plus = polynomial_continued_fraction_decomp(driving_point_Z_n_plus,driving_point_Z_d_plus)
    cfd_minus = polynomial_continued_fraction_decomp(driving_point_Z_n_minus,driving_point_Z_d_minus)

    s_coefs_plus, const_coef_plus = only_linear_term_of_cfd(cfd_plus)
    s_coefs_minus, const_coef_minus = only_linear_term_of_cfd(cfd_minus)

    Rload_plus = const_coef_plus
    Rload_minus = const_coef_minus

    C_plus = s_coefs_plus[1::2]
    L_plus = s_coefs_plus[::2]
    C_plus /= Rload_plus
    L_plus *= Rload_plus

    C_minus = s_coefs_minus[1::2]
    L_minus = s_coefs_minus[::2]
    C_minus /= Rload_minus
    L_minus *= Rload_minus

    ## Shunt first if L[0] < tol, so get rid of it
    shunt_first_plus = False
    if L_plus[0] < tol:
        L_plus = L_plus[1:]
        shunt_first_plus = True
    shunt_first_minus = False
    if L_minus[0] < tol:
        L_minus = L_minus[1:]
        shunt_first_minus = True

    if shunt_first_plus:
        lnf = LadderNetworkFilter(C_plus,L_plus,Rin=Rload_plus,Rout=Rload_plus,shunt_first=True)
        return lnf
    elif shunt_first_minus:
        lnf = LadderNetworkFilter(C_minus,L_minus,Rin=Rload_minus,Rout=Rload_minus,shunt_first=True)
        return lnf
    else:
        raise Exception(f"Both ladder networks are series first and that's not implemented")


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
    n = Polynomial(n).trim()
    d = Polynomial(d).trim()
    assert(len(n)<=len(d))

    # Only works for just 1 in the numerator for now
    assert(len(n) == 1)
    assert(n.coef[0] == 1.)

    range_mod_2 = np.arange(len(d)) % 2
    d_evens = Polynomial(d.coef * (range_mod_2 == 0)).trim()
    d_odds = Polynomial(d.coef * (range_mod_2 == 1)).trim()
    zfirst = len(d_evens) >= len(d_odds) # first component is impedance, otherwise is admittance (y)

    cfd = []
    if zfirst:
        cfd = polynomial_continued_fraction_decomp(d_evens,d_odds)
    else:
        cfd = polynomial_continued_fraction_decomp(d_odds,d_evens)

    for x in cfd:
        assert(len(x) == 2)
    cfd_s_coefs = np.array([x.coef[1] for x in cfd],dtype="float64")
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

    g = cauerI_synthesis_equal_inout_impedance([1],[1,2,2,1]) # 3rd order butter
    breakpoint()
    f = cauerI_synthesis_equal_inout_impedance([1],semi_gaussian_complex_all_pole_filter(3).to_tf().den)
    #h = cauerI_synthesis_equal_inout_impedance([1],[2,3,3,1])
    breakpoint()
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
