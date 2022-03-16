#!/usr/bin/env python3

import tempfile
from spiceAnalysis import SpiceAnalyzer
from spiceAnalysis import library as LIBRARY
import matplotlib.pyplot as mpl


class LadderNetworkFilter:
    """
    Class implementing a ladder network filter in SPICE
    """

    def __init__(self,Clist,Llist,Rin=50,Rout=50):
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

        self.order = len(Clist) + len(Llist)
        self.string = f".start {self.order}-Order Ladder Network Filter (autogen)\n"
        self.iFirstLCNode = 110
        self.iInputNode = 100
        self.iOutputNode = 199

        self.string += self._generate_ladder_string()
        self.string += self._generate_RinRout_string()
        self.string += ".end"

    def _generate_ladder_string(self):
        Clist = self.Clist
        Llist = self.Llist
        assert(len(Clist)>=len(Llist))
        assert(len(Clist)<=len(Llist)+1)
        assert(len(Clist) < self.iOutputNode-self.iFirstLCNode)

        result = ""

        for iC, C in enumerate(Clist):
            iNode = self.iFirstLCNode+iC
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
        with tempfile.TemporaryFile(mode="w+") as circuitFile:
            circuitFile.write(self.string)
            circuitFile.flush()
            circuitFile.seek(0)
            sa = SpiceAnalyzer(circuitFile)
            sa.analyzeFreqAndTrans(savename,title,self.iInputNode,0,self.iOutputNode,fstart,fstop,tstep,tstart,tstop,debug=debug)

if __name__ == "__main__":

    import numpy as np

    lnf = LadderNetworkFilter([1./50/2/np.pi,1./50/2/np.pi],[2*50./2/np.pi])
    #print(lnf.string)
    lnf.make_plots("LadderNetwork.pdf","3rd Order Butterworth Pi-Ladder Network",1e-3,1e3,1e-3,0,5,debug=False)
