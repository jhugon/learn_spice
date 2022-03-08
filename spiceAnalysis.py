#!/usr/bin/env python3

import re
import os
import tempfile
import subprocess
import io
import numpy
import matplotlib.pyplot as mpl

class TemplateModifier(object):

  def __init__(self,circuitTemplateFile):
    self.templateFile = circuitTemplateFile
    self.sourceLines = []
    self.optionLines = [".options noacct"]
    self.analyzerLines = []
    self.tmpFiles = []

  def addACSource(self,name,nodePlus,nodeMinus,mag,phase=0.):
    line = "{0} {1} {2} AC {3} {4}".format(name,nodePlus,nodeMinus,mag,phase)
    self.sourceLines.append(line)

  def addTransSource(self,name,nodePlus,nodeMinus,sourceStr):
    line = "{0} {1} {2} {3}".format(name,nodePlus,nodeMinus,sourceStr)
    self.sourceLines.append(line)

  def addACAnalysis(self,node,pointsPerDecade,fstart,fstop):
    line = ".ac dec {0:d} {1} {2}".format(pointsPerDecade,fstart,fstop)
    self.sourceLines.append(line)
    line = ".print ac {}".format(node)
    self.sourceLines.append(line)

  def addTransAnalysis(self,node,tstep,tstart,tstop):
    line = ".tran {0} {1} {2}".format(tstep,tstop,tstart)
    self.sourceLines.append(line)
    line = ".print tran {}".format(node)
    self.sourceLines.append(line)

  def clearSources(self):
    self.sourceLines.clear()

  def clearAnalyses(self):
    self.analyzerLines.clear()

  def clear(self):
    self.clearSources()
    self.clearAnalyses()
    self.tmpFiles.clear()

  def getFile(self):
    result = tempfile.NamedTemporaryFile(mode="w+",suffix=".cir",delete=True)
    if isinstance(self.templateFile,io.IOBase):
      self.templateFile.seek(0)
      for line in self.templateFile:
        if line[:4] != ".end" or line[:5] == ".ends":
          result.write(line)
    else:
      with open(self.templateFile) as infile:
          for line in infile:
            if line[:4] != ".end" or line[:5] == ".ends":
              result.write(line)
    for line in self.sourceLines:
      result.write(line + "\n")
    for line in self.optionLines:
      result.write(line + "\n")
    for line in self.analyzerLines:
      result.write(line + "\n")
    result.write(".end"+"\n")
    result.flush()
    self.tmpFiles.append(result)
    return result.name

  def __del__(self):
    for f in self.tmpFiles:
        f.close()

class SpiceAnalyzer(object):

  def __init__(self,circuitTemplateFile):
    self.circuitTemplateFile = circuitTemplateFile

  def getData(self,logfile):
    xtitle = None
    ytitles = None
    xdata = []
    ydata = []
    foundStart1 = False
    foundStart2 = False
    foundStart3 = False
    for line in logfile:
      if foundStart3:
        reStr = r"^(\d+)\s+([0-9eE.+-]+)\s+([0-9eE.+-]+)\s*$"
        datamatch = re.match(reStr,line)
        if datamatch:
          x = datamatch.group(2)
          x = float(x)
          xdata.append(x)
          y = datamatch.group(3)
          y = float(y)
          ydata.append(y)
      elif foundStart2:
        if re.match(r"^-+\s*$",line):
            foundStart3 = True
      elif foundStart1:
        reStr = r"^Index\s+(\w+)\s+([a-zA-Z0-9\(\)]+)\s*$"
        start2match = re.match(reStr,line)
        if start2match:
          xtitle = start2match.group(1)
          ytitle = start2match.group(2)
          foundStart2 = True
      elif re.match(r"^-+\s*$",line):
        foundStart1 = True
    xdata = numpy.array(xdata)
    ydata = numpy.array(ydata)
    return xdata, ydata, xtitle, ytitles


  def runSpice(self,circuitFileName,debug=False):
    if debug:
      print("runSpice: circuitFileName: '{}'".format(circuitFileName))
      with open(circuitFileName) as circuitFile:
          print(circuitFile.read())
    call_args = ["ngspice","-b",circuitFileName]
    output = tempfile.TemporaryFile(mode="w+")
    stderr = tempfile.TemporaryFile(mode="w+")
    try:
      subprocess.check_call(call_args,stdout=output,stderr=stderr)
    except subprocess.CalledProcessError as e:
      output.seek(0)
      stderr.seek(0)
      print("ngspice Error:")
      print(output.read())
      print(stderr.read())
    if debug:
      output.seek(0)
      print("output:")
      print(output.read())
    output.seek(0)
    xdata,ydata, xtitle, ytitle = self.getData(output)
    if debug:
      print(xdata)
      print(ydata)
    return xdata, ydata, xtitle, ytitle

  def analyzeAC(self,outFileName,inNodePlus,inNodeMinus,outProbes,mag,fstart,fstop,pointsPerDecade=5,current=False,debug=False):
    """
    outProbes is a list of things like '2' or '3,0' that can be put inside v() or vm().
    """
    xtitle = None
    freqGains = []
    magnitudes = []
    ytitles = []
    fig, axs = mpl.subplots(2)
    ax1 = axs[0]
    ax2 = axs[1]
    for outProbe in outProbes:
      outMagProbe = "vm({})".format(outProbe)
      template = TemplateModifier(self.circuitTemplateFile)
      template.addACSource("vac",inNodePlus,inNodeMinus,mag)
      template.addACAnalysis(outMagProbe,pointsPerDecade,fstart,fstop)
      circuitFileName = template.getFile()
      freqGain,magnitude, xtitle, ytitle = self.runSpice(circuitFileName,debug)
      freqGains.append(freqGain)
      magnitudes.append(magnitude)
      ytitles.append(ytitle)
    magnitudeDBs = [20*numpy.log10(ydata/mag) for ydata in magnitudes]
    freqPhases = []
    phases = []
    ytitles = []
    for outProbe in outProbes:
      outPhaseProbe = "vp({})".format(outProbe)
      template = TemplateModifier(self.circuitTemplateFile)
      if current:
        template.addACSource("iac",inNodePlus,inNodeMinus,mag)
      else:
        template.addACSource("vac",inNodePlus,inNodeMinus,mag)
      template.addACAnalysis(outPhaseProbe,pointsPerDecade,fstart,fstop)
      circuitFileName = template.getFile()
      xdata,ydata, xtitle, ytitle = self.runSpice(circuitFileName,debug)
      freqPhases.append(xdata)
      phases.append(ydata)
      ytitles.append(ytitle)
    phaseDegs = [ydata*180/numpy.pi for ydata in phases]
    for iCol in range(len(outProbes)):
      gainMargin, phaseMargin = self.getGainPhaseMargin(freqGains[iCol],magnitudeDBs[iCol],freqPhases[iCol],phaseDegs[iCol])
      label = "Point {}".format(outProbes[iCol])
      if not (gainMargin is None):
        label += ", GM: {:.0f} dB".format(gainMargin)
      if not (phaseMargin is None):
        label += r", PM: {:.0f}$^\circ$".format(phaseMargin)
      ax1.semilogx(freqGains[iCol],magnitudeDBs[iCol],label=label)
    #ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("V [dBc]")
    for iCol in range(len(outProbes)):
      ax2.semilogx(freqPhases[iCol],phaseDegs[iCol],label=label)
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel(r"Phase [$^\circ$]")
    ax1.legend(loc="best")
    if not(outFileName is None):
      fig.savefig(outFileName)
    return freqGains,magnitudeDBs,freqPhases,phaseDegs

  def getGainPhaseMargin(self,freqGain,magnitudeDB,freqPhase,phaseDeg):
    gainMargin = None
    phaseMargin = None
    unityFreq = self.findUnityGainFreq(freqGain,magnitudeDB)
    zeroPhaseFreq = self.findZeroPhaseFreq(freqPhase,phaseDeg)
    if unityFreq:
      phaseMargin = numpy.interp(unityFreq,freqPhase,phaseDeg)
    if zeroPhaseFreq:
      gainMargin = -numpy.interp(zeroPhaseFreq,freqGain,magnitudeDB)
    return gainMargin, phaseMargin

  def findUnityGainFreq(self,freqs,vDBc):
    """
    Finds the lowest frequency where the gain < 0 dBc.
    Returns None if it is the first freqency in the list.
    """
    if vDBc is None or len(vDBc) == 0:
      return None
    isBelowUnityGain = vDBc < 0
    result = None
    if not isBelowUnityGain[0]:
        freqsBelowUnityGain = freqs[isBelowUnityGain]
        if len(freqsBelowUnityGain) > 0:
            result = freqs[isBelowUnityGain][0]
    return result

  def findZeroPhaseFreq(self,freqs,phaseDeg):
    """
    Finds the lowest frequency where the phase changes sign.
    Returns None if one is not found
    """
    if phaseDeg is None or len(phaseDeg) == 0:
      return None
    signIsPos = phaseDeg > 0
    posFreqs = freqs[signIsPos]
    negFreqs = freqs[numpy.logical_not(signIsPos)]
    firstSignIsPos = signIsPos[0]
    result = None
    if firstSignIsPos:
      if len(negFreqs) > 0:
        result = negFreqs[0]
    else:
      if len(posFreqs) > 0:
        result = posFreqs[0]
    return result

  def analyzeTrans(self,outFileName,inNodePlus,inNodeMinus,outProbes,sourceStr,tstep,tstart,tstop,current=False,debug=False):
    """
    Transient Analysis
    outProbes is a list of things like '2' or '3,0' that can be put inside v() or vm().
    sourceStr is a str like PULSE(0,1,100ns,0,0,100ns,1s)
        where the pulse args are initial val, pulsed val, delay, rise time, fall time, pulse width, period.
    tstep is the recording step time
    tstart is the recording start time
    tstop is the recording stop time
    current: if true is a current source, else a voltage source
    """
    xtitle = None
    xdatas = []
    ydatas = []
    ytitles = []
    for outProbe in outProbes:
      outMagProbe = "v({})".format(outProbe)
      template = TemplateModifier(self.circuitTemplateFile)
      if current:
        template.addTransSource("itran",inNodePlus,inNodeMinus,sourceStr)
      else:
        template.addTransSource("vtran",inNodePlus,inNodeMinus,sourceStr)
      template.addTransAnalysis(outMagProbe,tstep,tstart,tstop)
      circuitFileName = template.getFile()
      xdata,ydata, xtitle, ytitle = self.runSpice(circuitFileName,debug)
      xdatas.append(xdata)
      ydatas.append(ydata)
      ytitles.append(ytitle)
    fig, ax = mpl.subplots()
    for iCol in range(len(outProbes)):
      label = outProbes[iCol]
      ax.plot(xdatas[iCol],ydatas[iCol],label=label)
    ax.set_xlabel(xtitle)
    ax.set_ylabel("V")
    ax.legend(loc="best")
    if not(outFileName is None):
      fig.savefig(outFileName)
    return xdatas,ydatas

  def analyzeManyTrans(self,outFileName,inNodePlus,inNodeMinus,outProbe,sourceStrs,tstep,tstart,tstop,current=False,debug=False):
    """
    Transient Analysis
    outProbes is a list of things like '2' or '3,0' that can be put inside v() or vm().
    sourceStrs is a list of strs like PULSE(0,1,100ns,0,0,100ns,1s)
        where the pulse args are initial val, pulsed val, delay, rise time, fall time, pulse width, period.
    tstep is the recording step time
    tstart is the recording start time
    tstop is the recording stop time
    current: if true is a current source, else a voltage source
    """
    xtitle = None
    xdatas = []
    ydatas = []
    ytitles = []
    for sourceStr in sourceStrs:
      outMagProbe = "v({})".format(outProbe)
      template = TemplateModifier(self.circuitTemplateFile)
      if current:
        template.addTransSource("itran",inNodePlus,inNodeMinus,sourceStr)
      else:
        template.addTransSource("vtran",inNodePlus,inNodeMinus,sourceStr)
      template.addTransAnalysis(outMagProbe,tstep,tstart,tstop)
      circuitFileName = template.getFile()
      xdata,ydata, xtitle, ytitle = self.runSpice(circuitFileName,debug)
      xdatas.append(xdata)
      ydatas.append(ydata)
      ytitles.append(ytitle)
    fig, ax = mpl.subplots()
    for iCol in range(len(sourceStrs)):
      label = sourceStrs[iCol]
      ax.plot(xdatas[iCol],ydatas[iCol],label=label)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(outProbe)
    ax.legend(loc="best")
    if not(outFileName is None):
      fig.savefig(outFileName)

def analyzeACManyCir(circuitTemplateFiles,outFileName,inNodePluses,inNodeMinuses,outProbes,labels,mag,fstart,fstop,pointsPerDecade=5,current=False,debug=False):
  freqGains = []
  magnitudeDBs = []
  freqPhases = []
  phaseDegs = []
  for circuitTemplateFile, inNodePlus, inNodeMinus, outProbe in zip(circuitTemplateFiles,inNodePluses, inNodeMinuses, outProbes):
    sa = SpiceAnalyzer(circuitTemplateFile)
    freqGain, magnitudeDB, freqPhase, phaseDeg = sa.analyzeAC(None,inNodePlus,inNodeMinus,outProbe,mag,fstart,fstop,pointsPerDecade=pointsPerDecade,current=current,debug=debug)
    if len(freqGains) > 1:
      raise Exception("freqGains should have length 1, freqGains: %s",freqGains)
    assert(len(freqGains) == len(freqPhases))
    freqGains.append(freqGain[0])
    magnitudeDBs.append(magnitudeDB[0])
    freqPhases.append(freqPhase[0])
    phaseDegs.append(phaseDeg[0])
  fig, axs = mpl.subplots(2)
  ax1 = axs[0]
  ax2 = axs[1]
  for iCol in range(len(circuitTemplateFiles)):
    #gainMargin, phaseMargin = self.getGainPhaseMargin(freqGains[iCol],magnitudeDBs[iCol],freqPhases[iCol],phaseDegs[iCol])
    label = labels[iCol]
    #if not (gainMargin is None):
    #  label += ", GM: {:.0f} dB".format(gainMargin)
    #if not (phaseMargin is None):
    #  label += r", PM: {:.0f}$^\circ$".format(phaseMargin)
    ax1.semilogx(freqGains[iCol],magnitudeDBs[iCol],label=label)
  #ax1.set_xlabel("Frequency [Hz]")
  ax1.set_ylabel("V [dBc]")
  for iCol in range(len(circuitTemplateFiles)):
    ax2.semilogx(freqPhases[iCol],phaseDegs[iCol],label=label)
  ax2.set_xlabel("Frequency [Hz]")
  ax2.set_ylabel(r"Phase [$^\circ$]")
  ax1.legend(loc="best")
  if not(outFileName is None):
    fig.savefig(outFileName)

def analyzeTransManyCir(circuitTemplateFiles,outFileName,inNodePluses,inNodeMinuses,outProbes,labels,sourceStr,tstep,tstart,tstop,current=False,debug=False):
  xdatas = []
  ydatas = []
  for circuitTemplateFile, inNodePlus, inNodeMinus, outProbe in zip(circuitTemplateFiles,inNodePluses, inNodeMinuses, outProbes):
    sa = SpiceAnalyzer(circuitTemplateFile)
    xdata, ydata = sa.analyzeTrans(None,inNodePlus,inNodeMinus,outProbe,sourceStr,tstep,tstart,tstop,current=current,debug=debug)
    if len(xdata) > 1:
      raise Exception("xdata should have length 1, xdata: %s",xdata)
    assert(len(xdata) == len(ydata))
    xdatas.append(xdata)
    ydatas.append(ydata)
  fig, ax = mpl.subplots()
  for iCol in range(len(circuitTemplateFiles)):
    print("lskdjf")
    ax.plot(xdatas[iCol],ydatas[iCol],label=labels[iCol])
  ax.set_xlabel("Time")
  ax.set_ylabel("V")
  #ax.legend(loc="best")
  if not(outFileName is None):
    fig.savefig(outFileName)

library = """*
*
* 1 input+, 2 input-, 3 output, 4 +V supply, 5 -V supply, 6 ground for all
*
* Ideal opamp
.subckt idealopamp 1 2 3 99 100 4
E1 3 4 1 2 1e8
.ends
*
* Ideal compensated opamp 1kHz GBW, 1e5 DC Gain
.subckt idealopamp1k 1 2 3 99 100 4
E5 5 4 1 2 1e5
r1 5 6 1e4
c1 6 4 1e-3
E6 3 4 6 4 1
.ends
*
* Ideal compensated opamp 10kHz GBW, 1e5 DC Gain
.subckt idealopamp10k 1 2 3 99 100 4
E5 5 4 1 2 1e5
r1 5 6 1e4
c1 6 4 1e-4
E6 3 4 6 4 1
.ends
*
* Ideal compensated opamp 100kHz GBW, 1e5 DC Gain
.subckt idealopamp100k 1 2 3 99 100 4
E5 5 4 1 2 1e5
r1 5 6 1e4
c1 6 4 1e-5
E6 3 4 6 4 1
.ends
*
* Ideal compensated opamp 1MHz GBW, 1e5 DC Gain
.subckt idealopamp1M 1 2 3 99 100 4
E5 5 4 1 2 1e5
r1 5 6 1e4
c1 6 4 1e-6
E6 3 4 6 4 1
.ends
*
* Ideal compensated opamp 10MHz GBW, 1e5 DC Gain
.subckt idealopamp10M 1 2 3 99 100 4
E5 5 4 1 2 1e5
r1 5 6 1e4
c1 6 4 1e-7
E6 3 4 6 4 1
.ends
*
* Ideal compensated opamp 100MHz GBW, 1e5 DC Gain
.subckt idealopamp100M 1 2 3 99 100 4
E5 5 4 1 2 1e5
r1 5 6 1e4
c1 6 4 1e-8
E6 3 4 6 4 1
.ends
*
* Ideal compensated opamp 1GHz GBW, 1e5 DC Gain
.subckt idealopamp1G 1 2 3 99 100 4
E5 5 4 1 2 1e5
r1 5 6 1e4
c1 6 4 1e-9
E6 3 4 6 4 1
.ends
*
* Ideal compensated opamp 10GHz GBW, 1e5 DC Gain
.subckt idealopamp10G 1 2 3 99 100 4
E5 5 4 1 2 1e5
r1 5 6 1e4
c1 6 4 1e-10
E6 3 4 6 4 1
.ends
*
*
"""

if __name__ == "__main__":
  sa = SpiceAnalyzer("example.cir.template")
  sa.analyzeAC("test.png",1,0,["2","3","4"],100,.01,10)
  sa.analyzeTrans("test2.png",1,0,["2","3","4"],"PULSE(0,1,0.1,0,0,0.1,100.)",0.01,0,3.)
  sa.analyzeManyTrans("test3.png",1,0,"2",
        [
            "PULSE(0,1,0.1,0,0,0.1,100.)",
            "PULSE(0,1,0.1,0,0,0.3,100.)",
            "PULSE(0,1,0.1,0,0,1.0,100.)",
            "PULSE(0,1,0.1,0,0,3.0,100.)",
        ],
        0.01,0,5.)

  with open("example.cir.template") as infile:
    sa = SpiceAnalyzer(infile)
    sa.analyzeAC("test4.png",1,0,["2","3","4"],100,.01,10)

  high = SpiceAnalyzer("high_pass.cir.template")
  high.analyzeAC("testHigh.png",1,0,["2"],1,"1k","1000k")
  high.analyzeManyTrans("testHigh2.png",1,0,"2",
        [
            #"PULSE(0,1,10u,0,0,1u,100u)",
            #"PULSE(0,1,10u,0,0,10u,100u)",
            "PULSE(0,1,10u,0,0,30u,100u)",
            "PULSE(0,1,10u,0,0,200u,100u)",
        ],
        "0.25u",0,"100u",debug=False)

  analyzeACManyCir(["example.cir.template","high_pass.cir.template"],"testMany.png",[1,1],[0,0],["4","2"],["example","high pass"],100,0.01,"1Meg",debug=False)
  analyzeTransManyCir(["example.cir.template","high_pass.cir.template"],"testManyTrans.png",[1,1],[0,0],["4","2"],["example","high pass"],"PULSE(0,1,0.1,0,0,3.0,100.)",0.01,0,5.,debug=False)
