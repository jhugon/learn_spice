#!/usr/bin/env python3

import tempfile
from spiceAnalysis import SpiceAnalyzer
from spiceAnalysis import library as LIBRARY
import matplotlib.pyplot as mpl
import numpy as np


class LadderNetworkFilter:
    """
    Class implementing a ladder network filter in SPICE
    """

    iFirstLCNode = 110
    iInputNode = 100
    iOutputNode = 199

    def __init__(self,Clist,Llist,Rin=50,Rout=50,shunt_first=True):
        """
        For now, only low pass is implemented and only shunt first

        Clist and Llist are lists of float capacitances and inductances, respectively, in F and H

        The input should be across 100 and 0 and the output accross 199 and 0

        A series input resistor and shunt output resistor are included
        """
        self.Clist = Clist
        self.Llist = Llist
        self.Rin = Rin
        self.Rout = Rout
        self.shunt_first=shunt_first

        self.order = len(Clist) + len(Llist)
        self.string = f".start {self.order}-Order Ladder Network Filter {'Shunt-first'*shunt_first}{'Series-first'*(not shunt_first)} (autogen)\n"

        self.string += self._generate_ladder_string()
        self.string += self._generate_RinRout_string()
        self.string += ".end"

        self.tempfile = None

    def __del__(self):
        try:
            self.tempfile.close()
        except AttributeError as e:
            pass
        except Exception as e:
            print(f"Exception closing tempfile: {type(e)} {e}")

    def __str__(self):
        result = "LadderNetworkFilter(C={},L={},Rin={},Rout={},shunt_first={})".format(self.Clist,self.Llist,self.Rin,self.Rout,self.shunt_first)
        return result

    def _generate_ladder_string(self):
        Clist = self.Clist
        Llist = self.Llist
        if self.shunt_first:
            assert(len(Clist)>=len(Llist))
            assert(len(Clist)<=len(Llist)+1)
            assert(len(Clist) < self.iOutputNode-self.iFirstLCNode)
        else:
            assert(len(Llist)>=len(Clist))
            assert(len(Llist)<=len(Clist)+1)
            assert(len(Llist) < self.iOutputNode-self.iFirstLCNode)

        result = ""

        for iC, C in enumerate(Clist):
            iNode = self.iFirstLCNode+iC
            if not self.shunt_first:
                iNode += 1
            result += f"C{iC} {iNode} 0 {C:g}\n"
        
        for iL, L in enumerate(Llist):
            iNode = self.iFirstLCNode+iL
            result += f"L{iL} {iNode} {iNode+1} {L:g}\n"
        return result

    def _generate_RinRout_string(self):
        result = f"Rin {self.iInputNode} {self.iFirstLCNode} {self.Rin}\n"
        result += f"Rout_wire {len(self.Llist)+self.iFirstLCNode} {self.iOutputNode} 1e-6\n"
        result += f"Rout {self.iOutputNode} 0 {self.Rout}\n"
        return result

    def make_plots(self,savename,title,fstart,fstop,tstep,tstart,tstop,debug=False):
        sa = self.get_spice_analyzer()
        sa.analyzeFreqAndTrans(savename,title,self.iInputNode,0,self.iOutputNode,fstart,fstop,tstep,tstart,tstop,debug=debug)

    def get_poles_zeros(self,debug=False):
        """
        returns (poles,zeros), where each is a numpy array of complex numbers
        """
        sa = self.get_spice_analyzer()
        return sa.analyzePolesZeros(None,self.iInputNode,0,self.iOutputNode,0,debug=debug)

    def get_spice_analyzer(self):
        try:
            self.tempfile.seek(0)
        except AttributeError as e:
            self.tempfile = tempfile.TemporaryFile(mode="w+")
            self.tempfile.write(self.string)
            self.tempfile.flush()
            self.tempfile.seek(0)
        return SpiceAnalyzer(self.tempfile)

    @classmethod
    def make_plots_many_filters(cls,filterList,labelList,savename,title,fstart,fstop,tstep,tstart,tstop,debug=False):
        SpiceAnalyzer.analyzeFreqAndTransManySpiceAnalyzers([x.get_spice_analyzer() for x in filterList],labelList,savename,title,cls.iInputNode,0,cls.iOutputNode,fstart,fstop,tstep,tstart,tstop,debug=debug)

def synchronouslyTunedFilter(n,f0,Q,R,shunt_first=True):
    """
    Makes a ladder filter with the given order, corner freq, Q, and Rin/out
    """

    assert(n>1)
    factor = np.sqrt(2**(1/n)-1)
    fsection = f0*factor**(-2)
    Qsection = Q*factor
    L = R/(2*np.pi*Qsection*fsection)
    C = Qsection/(2*np.pi*R*fsection)
    nC = n//2 + n % 2
    nL = n//2
    result = LadderNetworkFilter([C]*nC,[L]*nL,Rin=R,Rout=R,shunt_first=shunt_first)
    return result

if __name__ == "__main__":


    butterworth3 = LadderNetworkFilter([1./50/2/np.pi,1./50/2/np.pi],[2*50./2/np.pi])
    #butterworth3.make_plots("LadderNetwork.pdf","3rd Order Butterworth Pi-Ladder Network",1e-3,1e3,1e-3,0,5,debug=False)
    #butterworth3.get_spice_analyzer().analyzeZinZout("ZinZout.pdf",100,0,199,0,1e-3,1e3)

    bessel2 = LadderNetworkFilter([0.5755/50/2/np.pi],[2.1478*50./2/np.pi])
    bessel3 = LadderNetworkFilter([0.3374/50/2/np.pi,2.2034/50/2/np.pi],[0.9705*50./2/np.pi])
    bessel5 = LadderNetworkFilter([0.1743/50/2/np.pi,0.8040/50/2/np.pi,2.2582/50/2/np.pi],[0.5072*50./2/np.pi,1.1110*50./2/np.pi])
    
    # for a single L and single C, C = 1/(4*pi*R*f_0) and L = R/(pi*f_0) is critically damped
    TotalC = 0.5
    TotalL = 2
    TotalC /= 50*2*np.pi
    TotalL /= 2*np.pi/50
    simpleCLR = LadderNetworkFilter([TotalC],[TotalL])
    simpleLCR = LadderNetworkFilter([TotalC],[TotalL],shunt_first=False)
    simpleT = LadderNetworkFilter([TotalC],[TotalL/2,TotalL/2],shunt_first=False)
    simplePi = LadderNetworkFilter([TotalC*2,TotalC*2],[TotalL])

    R = 50.
    Q = 0.5
    f0 = 1
    synchPi2 = synchronouslyTunedFilter(2,f0,Q,R)
    synchPi3 = synchronouslyTunedFilter(3,f0,Q,R)
    synchPi4 = synchronouslyTunedFilter(4,f0,Q,R)
    synchPi5 = synchronouslyTunedFilter(5,f0,Q,R)

    poles, zeros = synchPi3.get_poles_zeros()
    print("poles:",poles)

    filtersAndLabels = [
        #(butterworth3,"Butterworth 3O"),
        #(bessel2,"Bessel 2O"),
        #(bessel3,"Bessel 3O"),
        #(bessel5,"Bessel 5O"),
        #(simpleCLR,"Simple CLR"),
        #(simpleLCR,"Simple LCR"),
        #(simpleT,"Simple T"),
        #(simplePi,"Simple $\Pi$"),
        (synchPi2,"Synch $\Pi$ 2O"),
        (synchPi3,"Synch $\Pi$ 3O"),
        #(synchPi4,"Synch $\Pi$ 4O"),
        #(synchPi5,"Synch $\Pi$ 5O"),
    ]
    filters = [x[0] for x in filtersAndLabels]
    filterLabels = [x[1] for x in filtersAndLabels]
    LadderNetworkFilter.make_plots_many_filters(filters,filterLabels,"LadderFilters.pdf","Comparison of Ladder Filters",1e-3,1e3,1e-3,0,5,debug=False)

    trySynthPi = LadderNetworkFilter([2/3.,1/3.],[9/2.],Rin=1.,Rout=1.)
    trySynthPi.make_plots("TrySynthPi.pdf",r"Sythesized $\Pi$ filter for $\frac{1}{(s+1)^3}$",1e-4,1e3,13-3,0,100)
    print(trySynthPi.string)
    print(trySynthPi.get_poles_zeros())
    trySynthPi2 = LadderNetworkFilter([2.],[0.5],Rin=1.,Rout=1.)
    trySynthPi2.make_plots("TrySynthPi2.pdf",r"Sythesized $\Pi$ filter for $\frac{1}{(s+1)^2}$",1e-4,1e3,13-3,0,100)
    print(trySynthPi2.get_poles_zeros())
