#!/usr/bin/env python3

from scipy import signal, special, cluster, optimize
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np

def plot_filters_behavior(systems,labels,title,savename,t_max=7,f_min=1e-3,f_max=1e3):
    N_ts = 200
    N_ws = 200
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

def semi_gaussian_complex_pole_locations(N,out_img_fn=None):
    """
    The poles are the complex roots of the equation:

    0 = 1 - Sum_1^N s^{2*i}/i!
    """

    def eqn(s):
        result = 1.+0j
        for i in range(1,N+1):
            result -= s**i/special.factorial(i)
        return result

    a = np.linspace(-5.5,2.5,300)
    b = np.linspace(-5.5,5.5,300)
    av,bv = np.meshgrid(a,b)
    cv = av+bv*1j
    cf = eqn(cv)
    mf = np.abs(cf)

    agood = av[mf < 0.1]
    bgood = bv[mf < 0.1]
    obsgood = np.column_stack([agood,bgood])
    centroids2, labels = cluster.vq.kmeans2(obsgood,N,minit="++")

    if out_img_fn:
        fig, ax = plt.subplots(constrained_layout=True)
        pcm = ax.pcolormesh(av,bv,mf,shading="auto",norm=mcolors.LogNorm(vmax=10))
        #ax.scatter(agood,bgood,marker='x',c='r')
        ax.scatter(centroids2[:,0],centroids2[:,1],marker='x',c='m')
        ax.set_xlabel("Re(s)")
        ax.set_ylabel("Im(s)")
        fig.colorbar(pcm)
        fig.savefig("GaussianPoleLocations.png")

    poles = np.zeros(N,dtype="complex64")
    for iPole in range(centroids2.shape[0]):
        centroid = centroids2[iPole]
        x0 = centroid[0]+centroid[1]*1j
        x1 = x0 + 0.05+0.05j
        solution = optimize.root_scalar(eqn,x0=x0,x1=x1)
        if not solution.converged:
            raise Exception(f"Root finding did not converge {solution}")
        poles[iPole] = solution.root
    return poles
    

def semi_gaussian_complex_all_pole_filter(N):
    poles = semi_gaussian_complex_pole_locations(N)
    poles *= 0.2
    f = signal.ZerosPolesGain([],poles,[1])
    _, resp = f.freqresp([1e-9])
    scale_factor = 1./abs(resp)
    f = signal.ZerosPolesGain([],poles,scale_factor)
    return f


if __name__ == "__main__":
    alpha = -1
    n = 2
    ## The peak of a real all-same-real-pole filter is delayed by (n-1)/(-alpha)
    ## Critically damped is all real!
    normalized_alpha = 1-n
    #real_all_pole_filter = signal.ZerosPolesGain([],[alpha]*n,[1])
    real_all_pole_filter = signal.ZerosPolesGain([],[normalized_alpha]*n,[1])

    complex_pole_filter = signal.ZerosPolesGain([],[-2+0.4j]*3,[1])

    # From https://www.bnl.gov/tcp/uploads/files/BSA11-10j.pdf "Shaper Design
    # in CMOS for High Dynamic Range" by G De Geronimo and S Li
    shaping_time_sf = 2
    p0 = -1.793/shaping_time_sf
    w1 = 1.976/shaping_time_sf
    Q1 = 0.606
    zeta1 = 0.5/Q1
    p1 =-zeta1*w1 + 1j*w1*np.sqrt(1-zeta1**2)
    semi_gaussian_C3_filter = signal.ZerosPolesGain([],[p0,p1],[1])
    print(semi_gaussian_C3_filter.to_tf())

    filters_to_plot = [
        (real_all_pole_filter,"Real All-Pole Filter"),
        (complex_pole_filter,"Complex Pole Filter"),
        (semi_gaussian_C3_filter,"Semi-Gaussian C3 Filter (from paper)"),
        (semi_gaussian_complex_all_pole_filter(2),"Semi-Gaussian C2 Filter"),
        #(semi_gaussian_complex_all_pole_filter(3),"Semi-Gaussian C3 Filter"),
        #(semi_gaussian_complex_all_pole_filter(4),"Semi-Gaussian C4 Filter"),
        #(semi_gaussian_complex_all_pole_filter(5),"Semi-Gaussian C5 Filter"),
    ]
    plot_filters_behavior([x[0] for x in filters_to_plot],[x[1] for x in filters_to_plot],"Filter Design","Filter_Design.pdf",t_max=15)
