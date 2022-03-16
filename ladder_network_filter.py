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
        impulseString = f"PULSE(0,1,0,{tstep},{tstep},{tstep},{tstop*1000})"
        stepString = f"PWL(0,0,{tstep},1)"
        with tempfile.TemporaryFile(mode="w+") as circuitFile:
            circuitFile.write(self.string)
            circuitFile.flush()
            circuitFile.seek(0)
            sa = SpiceAnalyzer(circuitFile)
            freqGains,magnitudeDBs,freqPhases,phaseDegs = sa.analyzeAC(None,self.iInputNode,0,[self.iOutputNode],1,fstart,fstop,debug=debug)
            tImpulse, vImpulse = sa.analyzeTrans(None,self.iInputNode,0,[self.iOutputNode],impulseString,tstep,tstart,tstop,debug=debug)
            tStep, vStep = sa.analyzeTrans(None,self.iInputNode,0,[self.iOutputNode],stepString,tstep,tstart,tstop,debug=debug)
            fig, ((ax_g,ax_i),(ax_p,ax_s)) = mpl.subplots(nrows=2,ncols=2,figsize=(8.5,11),constrained_layout=True,sharex="col")
            ax_g.plot(freqGains[0],magnitudeDBs[0],label="Amplitude")
            ax_p.plot(freqPhases[0],phaseDegs[0],label="Phase")
            ax_p.set_xlabel("Frequency [Hz]")
            ax_g.set_ylabel("Voltage Gain/Loss [dB]")
            ax_p.set_ylabel("Phase [deg]")
            ax_g.set_xscale("log")
            ax_p.set_xscale("log")
            ax_g.set_xlim(fstart,fstop)
            ax_p.set_xlim(fstart,fstop)
            ax_g.tick_params(axis='both',which='both',direction="in",bottom=True,top=True,left=True,right=True)
            ax_p.tick_params(axis='both',which='both',direction="in",bottom=True,top=True,left=True,right=True)
            ax_i.plot(tImpulse[0],vImpulse[0],label="Impulse")
            ax_s.plot(tStep[0],vStep[0],label="Step")
            ax_s.set_xlabel("t [s]")
            ax_i.set_ylabel("Impulse Response [V/V]")
            ax_s.set_ylabel("Step Response [V/V]")
            ax_i.set_xlim(tstart,tstop)
            ax_s.set_xlim(tstart,tstop)
            ax_i.tick_params(axis='both',which='both',direction="in",bottom=True,top=True,left=True,right=True)
            ax_s.tick_params(axis='both',which='both',direction="in",bottom=True,top=True,left=True,right=True)
            fig.suptitle(title)
            fig.savefig(savename)

if __name__ == "__main__":

    lnf = LadderNetworkFilter([1,1],[2])
    print(lnf.string)
    lnf.make_plots("LadderNetwork.pdf","3rd Order Ladder Network",1e-3,1e3,1e-3,0,20,debug=False)
