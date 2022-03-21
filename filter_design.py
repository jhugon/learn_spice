#!/usr/bin/env python3

from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

def plot_filter_behavior(system,savename):
    N_ts = 200
    t_max = 7
    ts = np.linspace(0,t_max,N_ts)
    t_impulse, y_impulse = signal.impulse(system,T=ts)
    t_step, y_step = signal.step(system,T=ts)

    w, mag, phase = signal.bode(system)

    fig, ((ax_f, ax_i), (ax_p, ax_s)) = plt.subplots(nrows=2,ncols=2,figsize=(11.,8.5),constrained_layout=True,sharex="col")

    ax_i.plot(t_impulse,y_impulse,label="Impulse")
    ax_s.plot(t_step,y_step,label="Step")
    ax_s.set_xlabel("t [rad?]")
    ax_i.set_ylabel("Impulse Response [V/V]")
    ax_s.set_ylabel("Step Response [V/V]")
    #for i in range(10):
    #    ax_t.axvline(i,ls='--',c="0.5")
    ax_i.axhline(0.,ls='--',c="0.5")
    ax_i.set_xlim(0,t_max)
    ax_s.set_xlim(0,t_max)

    ax_f.plot(w,mag,label="Magnitude")
    ax_p.plot(w,phase,label="Phase")
    ax_p.set_xlabel("f [rad/s]")
    ax_f.set_ylabel("Magnitude [dB]")
    ax_p.set_ylabel("Phase [deg]")
    #ax_f.set_ylim(-180,0)
    ax_f.axhline(0.,ls='--',c="0.5")
    ## Make phase labels at multiples of 45 deg
    ax_p_yticks = np.arange(phase.min()//45,np.ceil(phase.max()/45.)+1)*45
    ax_p.set_ylim(ax_p_yticks[0],ax_p_yticks[-1])
    ax_p.set_yticks(ax_p_yticks)

    ax_p.set_xlim(w[0],w[-1])
    ax_f.set_xlim(w[0],w[-1])


    fig.savefig(savename)

if __name__ == "__main__":
    alpha = -1
    n = 2
    ## The peak of a real all-same-real-pole filter is delayed by (n-1)/(-alpha)
    ## Critically damped is all real!
    normalized_alpha = 1-n
    #real_all_pole_filter = signal.ZerosPolesGain([],[alpha]*n,[1])
    real_all_pole_filter = signal.ZerosPolesGain([],[normalized_alpha]*n,[1])
    plot_filter_behavior(real_all_pole_filter,"real_all_pole_filter.pdf")
    import sys
    sys.exit(0)

    complex_pole_filter = signal.ZerosPolesGain([],[-2+0.4j]*3,[1])
    plot_filter_behavior(complex_pole_filter,"complex_pole_filter.pdf")

    # From https://www.bnl.gov/tcp/uploads/files/BSA11-10j.pdf "Shaper Design
    # in CMOS for High Dynamic Range" by G De Geronimo and S Li
    shaping_time_sf = 2
    p0 = -1.793/shaping_time_sf
    w1 = 1.976/shaping_time_sf
    Q1 = 0.606
    zeta1 = 0.5/Q1
    p1 =-zeta1*w1 + 1j*w1*np.sqrt(1-zeta1**2)
    semi_gaussian_C3_filter = signal.ZerosPolesGain([],[p0,p1],[1])
    plot_filter_behavior(semi_gaussian_C3_filter,"semi_gaussian_C3_filter.pdf")
    print(semi_gaussian_C3_filter.to_tf())

