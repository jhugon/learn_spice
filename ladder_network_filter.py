#!/usr/bin/env python3

import tempfile
from spiceAnalysis import SpiceAnalyzer
from spiceAnalysis import library as LIBRARY
import matplotlib.pyplot as mpl


class LadderNetworkFilter:
    """
    Class implementing a ladder network filter in SPICE
    """

    iFirstLCNode = 110
    iInputNode = 100
    iOutputNode = 199

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
        sa = self.get_spice_analyzer()
        sa.analyzeFreqAndTrans(savename,title,self.iInputNode,0,self.iOutputNode,fstart,fstop,tstep,tstart,tstop,debug=debug)

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


if __name__ == "__main__":

    import numpy as np

    butterworth3 = LadderNetworkFilter([1./50/2/np.pi,1./50/2/np.pi],[2*50./2/np.pi])
    butterworth3.make_plots("LadderNetwork.pdf","3rd Order Butterworth Pi-Ladder Network",1e-3,1e3,1e-3,0,5,debug=False)
    butterworth3.get_spice_analyzer().analyzeZinZout("ZinZout.pdf",100,0,199,0,1e-3,1e3)

    bessel3 = LadderNetworkFilter([0.3374/50/2/np.pi,2.2034/50/2/np.pi],[0.9705*50./2/np.pi])
    bessel5 = LadderNetworkFilter([0.1743/50/2/np.pi,0.8040/50/2/np.pi,2.2582/50/2/np.pi],[0.5072*50./2/np.pi,1.1110*50./2/np.pi])

    filters = [butterworth3,bessel3,bessel5]
    filterLabels = ["Butterworth 3O","Bessel 3O","Bessel 5O"]
    LadderNetworkFilter.make_plots_many_filters(filters,filterLabels,"LadderFilters.pdf","Comparison of Ladder Filters",1e-3,1e3,1e-3,0,3,debug=False)
