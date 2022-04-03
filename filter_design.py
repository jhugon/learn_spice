#!/usr/bin/env python3

from scipy import signal, special, cluster, optimize
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np

def plot_filters_behavior(systems,labels,title,savename,t_max=7,f_min=1e-3,f_max=1e3):
    N_ts = 500
    N_ws = 500
    ws = np.logspace(np.log10(f_min),np.log10(f_max),N_ws)
    ts = np.linspace(0,t_max,N_ts)

    system_results = []
    for system,label in zip(systems,labels):
        _, y_impulse = signal.impulse(system,T=ts)
        _, y_step = signal.step(system,T=ts)
        w, mag, phase = signal.bode(system,w=ws)
        system_results.append((label,y_impulse,y_step,w,mag,phase))

    fig, ((ax_f, ax_i), (ax_p, ax_s)) = plt.subplots(nrows=2,ncols=2,figsize=(11.,8.5),constrained_layout=True,sharex="col")
    ax_f.set_xscale("log")
    ax_p.set_xscale("log")
    ax_p.set_xlim(ws[0],ws[-1])
    ax_f.set_xlim(ws[0],ws[-1])
    ax_s.set_xlabel("t [rad?]")
    ax_i.set_ylabel("Impulse Response [V/V]")
    ax_s.set_ylabel("Step Response [V/V]")
    #for i in range(10):
    #    ax_t.axvline(i,ls='--',c="0.5")
    ax_i.axhline(0.,ls='--',c="0.5")
    ax_i.set_xlim(0,t_max)
    ax_s.set_xlim(0,t_max)
    ax_p.set_xlabel("f [rad/s]")
    ax_f.set_ylabel("Magnitude [dB]")
    ax_p.set_ylabel("Phase [deg]")
    #ax_f.set_ylim(-180,0)
    #ax_f.axhline(0.,ls='--',c="0.5")


    for label,y_impulse,y_step,w,mag,phase in system_results:
        ax_i.plot(ts,y_impulse,label=label)
        ax_s.plot(ts,y_step,label=label)
        ax_f.plot(ws,mag,label=label)
        ax_p.plot(ws,phase,label=label)

    ## Make phase labels at multiples of 45 deg
    ax_p_yticks = np.arange(min([x[5].min() for x in system_results])//45,np.ceil(max([x[5].max() for x in system_results])/45.)+1)*45
    ax_p.set_ylim(ax_p_yticks[0],ax_p_yticks[-1])
    ax_p.set_yticks(ax_p_yticks)

    if len(system_results) > 1:
        ax_f.legend()

    fig.suptitle(title)

    fig.savefig(savename)

    print("Frequency Response:")
    for label,y_impulse,y_step,w,mag,phase in system_results:
        w_3db = float('nan')
        f_3db = float('nan')
        w_m180 = float('nan')
        f_m180 = float('nan')
        mag_phase_m180 = float('nan')
        gain_margin = float('nan')
        try:
            w_3db = ws[mag <= mag[0]-3.][0]
            f_3db = w_3db/2/np.pi
        except IndexError:
            pass
        try:
            w_m180 = ws[phase <= -180.][0]
            f_m180 = w_m180/2/np.pi
        except IndexError:
            pass
        try:
            mag_phase_m180 = mag[phase <= -180.][0]
            gain_margin = -mag_phase_m180
        except IndexError:
            pass
        print(f"{label:50} -3db: {w_3db:8.3g} rad/s = {f_3db:8.3g} Hz,  phase<=180 @ {w_m180:8.3g} rad/s = {f_m180:8.3g} Hz,  gain margin: {gain_margin:4.1f} dB")
    print("Impulse Response:")
    for label,y_impulse,y_step,w,mag,phase in system_results:
        max_impulse = float('nan')
        min_impulse = float('nan')
        iMax_impulse = float('nan')
        tMax_impulse = float('nan')
        t10pct_impulse = float('nan')
        t100pct_impulse = float('nan')
        t10pct_fall_impulse = float('nan')
        fwhm_impulse = float('nan')
        try:
            max_impulse = y_impulse.max()
            min_impulse = y_impulse.min()
            iMax_impulse = y_impulse.argmax()
        except ValueError:
            pass
        try:
            tMax_impulse = ts[iMax_impulse]
        except IndexError:
            pass
        try:
            t10pct_impulse = ts[:iMax_impulse][y_impulse[:iMax_impulse] <= max_impulse*0.1][-1]
        except IndexError:
            pass
        try:
            t100pct_impulse = ts[y_impulse == max_impulse][0]
        except IndexError:
            pass
        try:
            t10pct_fall_impulse = ts[iMax_impulse:][y_impulse[iMax_impulse:] <= max_impulse*0.1][0]
        except IndexError:
            pass
        try:
            rise_half_impulse = ts[:iMax_impulse][y_impulse[:iMax_impulse] <= max_impulse*0.5][-1]
            fall_half_impulse = ts[iMax_impulse:][y_impulse[iMax_impulse:] <= max_impulse*0.5][0]
            fwhm_impulse = fall_half_impulse-rise_half_impulse
        except IndexError:
            pass
        print(f"{label:50} peak delay: {tMax_impulse:8.3g}  10-100%: {t100pct_impulse-t10pct_impulse:8.3g}  100-10%: {t10pct_fall_impulse-t100pct_impulse:5.3g}  max: {max_impulse:5.3g}  min: {min_impulse:4.3g}  FWHM: {fwhm_impulse:8.5g}")
    print("Step Response:")
    for label,y_impulse,y_step,w,mag,phase in system_results:
        max_step = float('nan')
        overshoot = float('nan')
        t10pct_step = float('nan')
        t90pct_step = float('nan')
        try:
            max_step = y_step.max()
        except ValueError:
            pass
        try:
            overshoot = max_step - y_step[-1]
        except IndexError:
            pass
        try:
            t10pct_step = ts[y_step >= 0.1*y_step[-1]][0]
        except IndexError:
            pass
        try:
            t90pct_step = ts[y_step >= 0.9*y_step[-1]][0]
        except IndexError:
            pass
        print(f"{label:50} 10-90%: {t90pct_step-t10pct_step:8.4g}  overshoot: {overshoot:8.4g}")

def semi_gaussian_complex_pole_locations(N,out_img_fn=None):
    """
    The poles are the complex roots of the equation:

    0 = Sum_1^N (-1)^i*s^{2*i}/i!
    """

    assert(N >= 2)

    poly_coefs = np.zeros(2*N+1)
    for j in range(2*N+1):
        i = j // 2
        if j % 2 != 0:
            continue
        poly_coefs[j] = (-1.)**i/special.factorial(i)
    poly = np.polynomial.Polynomial(poly_coefs)
    poly_roots = poly.roots()
    poles = poly_roots[poly_roots.real <= 0.]
    # normalize to a peak delay of 1 s
    #if N == 2:
    #    poles *= 0.8657314629
    #elif N == 3:
    #    poles *= 1.482965932
    #elif N == 4:
    #    poles *= 1.955911824
    #elif N == 5:
    #   poles *= 2.340681363
    return poles
    

def semi_gaussian_complex_all_pole_filter(N,out_img_fn=None):
    poles = semi_gaussian_complex_pole_locations(N)
    f1 = signal.ZerosPolesGain([],poles,[1])

    return f1
#    ts, y_impulse = f1.impulse()
#    iMax_impulse = y_impulse.argmax()
#    poles_scale_factor = ts[iMax_impulse]
#    poles *= poles_scale_factor
#
#    f2 = signal.ZerosPolesGain([],poles,[1])
#    _, resp = f2.freqresp([1e-9])
#    scale_factor = 1./abs(resp)
#
#    f3 = signal.ZerosPolesGain([],poles,scale_factor)
#    return f3


if __name__ == "__main__":

    alpha = -1
    n = 2
    ## The peak of a real all-same-real-pole filter is delayed by (n-1)/(-alpha)
    ## Critically damped is all real!
    normalized_alpha = 1-n
    #real_all_pole_filter = signal.ZerosPolesGain([],[alpha]*n,[1])
    real_all_pole_filter = signal.ZerosPolesGain([],[normalized_alpha]*n,[1])

    complex_pole_filter = signal.ZerosPolesGain([],[-2+0.4j]*3,[8.48476847])

    # From https://www.bnl.gov/tcp/uploads/files/BSA11-10j.pdf "Shaper Design
    # in CMOS for High Dynamic Range" by G De Geronimo and S Li
    shaping_time_sf = 2
    p0 = -1.793/shaping_time_sf
    w1 = 1.976/shaping_time_sf
    Q1 = 0.606
    zeta1 = 0.5/Q1
    p1 =-zeta1*w1 + 1j*w1*np.sqrt(1-zeta1**2)
    semi_gaussian_C3_filter = signal.ZerosPolesGain([],[p0,p1],[1])

    for i in range(2,12):
        semi_gaussian_complex_pole_locations(i,f"GaussianPoleLocations_{i}.png")

    filters_to_plot = [
        #(real_all_pole_filter,"Real All-Pole Filter"),
        #(complex_pole_filter,"Complex Pole Filter"),
        #(semi_gaussian_C3_filter,"Semi-Gaussian C3 Filter (from paper)"),
        #(semi_gaussian_complex_all_pole_filter(2),"Semi-Gaussian C2 Filter"),
        #(semi_gaussian_complex_all_pole_filter(3),"Semi-Gaussian C3 Filter"),
        #(semi_gaussian_complex_all_pole_filter(4),"Semi-Gaussian C4 Filter"),
        #(semi_gaussian_complex_all_pole_filter(5),"Semi-Gaussian C5 Filter"),
        #(semi_gaussian_complex_all_pole_filter(6),"Semi-Gaussian C6 Filter"),
        #(semi_gaussian_complex_all_pole_filter(7),"Semi-Gaussian C7 Filter"),
        #(semi_gaussian_complex_all_pole_filter(8),"Semi-Gaussian C8 Filter"),
        #(semi_gaussian_complex_all_pole_filter(9),"Semi-Gaussian C9 Filter"),
        #(semi_gaussian_complex_all_pole_filter(10),"Semi-Gaussian C10 Filter"),
        #(semi_gaussian_complex_all_pole_filter(10),"Semi-Gaussian C10 Filter"),
        (signal.TransferFunction(*signal.bessel(2,1,btype="lowpass",analog=True,output="ba")),"Bessel 2O Filter"),
        (signal.TransferFunction(*signal.bessel(3,1,btype="lowpass",analog=True,output="ba")),"Bessel 3O Filter"),
        (signal.TransferFunction(*signal.bessel(4,1,btype="lowpass",analog=True,output="ba")),"Bessel 4O Filter"),
        (signal.TransferFunction(*signal.bessel(5,1,btype="lowpass",analog=True,output="ba")),"Bessel 5O Filter"),
    ]
    plot_filters_behavior([x[0] for x in filters_to_plot],[x[1] for x in filters_to_plot],"Filter Design","Filter_Design.pdf",t_max=15)
