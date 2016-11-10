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
    line = "{0} {1:d} {2:d} AC {3} {4}".format(name,nodePlus,nodeMinus,mag,phase)
    self.sourceLines.append(line)

  def addTransSource(self,name,nodePlus,nodeMinus,sourceStr):
    line = "{0} {1:d} {2:d} {3}".format(name,nodePlus,nodeMinus,sourceStr)
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
    result = tempfile.NamedTemporaryFile(mode="w+")
    if isinstance(self.templateFile,io.IOBase):
      self.templateFile.seek(0)
      for line in self.templateFile:
        if line[:4] != ".end":
          result.write(line)
    else:
      with open(self.templateFile) as infile:
          for line in infile:
            if line[:4] != ".end":
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
      print(output.read())
    output.seek(0)
    xdata,ydata, xtitle, ytitle = self.getData(output)
    if debug:
      print(xdata)
      print(ydata)
    return xdata, ydata, xtitle, ytitle

  def analyzeAC(self,outFileName,inNodePlus,inNodeMinus,outProbes,mag,fstart,fstop,pointsPerDecade=5,debug=False):
    """
    outProbes is a list of things like 'v(2)', 'vm(3,0)'
    """
    xdata = None
    xtitle = None
    ydatas = []
    ytitles = []
    for outProbe in outProbes:
      template = TemplateModifier(self.circuitTemplateFile)
      template.addACSource("vac",inNodePlus,inNodeMinus,mag)
      template.addACAnalysis(outProbe,pointsPerDecade,fstart,fstop)
      circuitFileName = template.getFile()
      xdata,ydata, xtitle, ytitle = self.runSpice(circuitFileName,debug)
      ydatas.append(ydata)
      ytitles.append(ytitle)
    ydatas = tuple(ydatas)
    ydata = numpy.stack(tuple(ydatas))
    ydataDB = 20*numpy.log(ydata/mag)
    fig, ax = mpl.subplots()
    for iCol in range(len(outProbes)):
      label = outProbes[iCol]
      ax.semilogx(xdata,ydataDB[iCol],label=label)
    ax.set_xlabel(xtitle)
    ax.set_ylabel("V [dBc]")
    ax.legend()
    fig.savefig(outFileName)

  def analyzeTrans(self,outFileName,inNodePlus,inNodeMinus,outProbes,sourceStr,tstep,tstart,tstop,debug=False):
    """
    Transient Analysis
    outProbes is a list of things like 'v(2)', 'vm(3,0)'
    sourceStr is a str like PULSE(0,1,100ns,0,0,100ns,1s)
        where the pulse args are initial val, pulsed val, delay, rise time, fall time, pulse width, period.
    tstep is the recording step time
    tstart is the recording start time
    tstop is the recording stop time
    """
    xdata = None
    xtitle = None
    ydatas = []
    ytitles = []
    for outProbe in outProbes:
      template = TemplateModifier(self.circuitTemplateFile)
      template.addTransSource("vtran",inNodePlus,inNodeMinus,sourceStr)
      template.addTransAnalysis(outProbe,tstep,tstart,tstop)
      circuitFileName = template.getFile()
      xdata,ydata, xtitle, ytitle = self.runSpice(circuitFileName,debug)
      ydatas.append(ydata)
      ytitles.append(ytitle)
    ydatas = tuple(ydatas)
    ydata = numpy.stack(tuple(ydatas))
    fig, ax = mpl.subplots()
    for iCol in range(len(outProbes)):
      label = outProbes[iCol]
      ax.plot(xdata,ydata[iCol],label=label)
    ax.set_xlabel(xtitle)
    ax.set_ylabel("V")
    ax.legend()
    fig.savefig(outFileName)

  def analyzeManyTrans(self,outFileName,inNodePlus,inNodeMinus,outProbe,sourceStrs,tstep,tstart,tstop,debug=False):
    """
    Transient Analysis
    outProbe is a str like 'v(2)', 'vm(3,0)'
    sourceStrs is a list of strs like PULSE(0,1,100ns,0,0,100ns,1s)
        where the pulse args are initial val, pulsed val, delay, rise time, fall time, pulse width, period.
    tstep is the recording step time
    tstart is the recording start time
    tstop is the recording stop time
    """
    xdata = None
    xtitle = None
    ydatas = []
    ytitles = []
    for sourceStr in sourceStrs:
      template = TemplateModifier(self.circuitTemplateFile)
      template.addTransSource("vtran",inNodePlus,inNodeMinus,sourceStr)
      template.addTransAnalysis(outProbe,tstep,tstart,tstop)
      circuitFileName = template.getFile()
      xdata,ydata, xtitle, ytitle = self.runSpice(circuitFileName,debug)
      ydatas.append(ydata)
      ytitles.append(ytitle)
    ydatas = tuple(ydatas)
    ydata = numpy.stack(tuple(ydatas))
    fig, ax = mpl.subplots()
    for iCol in range(len(sourceStrs)):
      label = sourceStrs[iCol]
      ax.plot(xdata,ydata[iCol],label=label)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(outProbe)
    ax.legend()
    fig.savefig(outFileName)

if __name__ == "__main__":
  sa = SpiceAnalyzer("example.cir.template")
  sa.analyzeAC("test.png",1,0,["vm(2)","vm(3)","vm(4)"],100,.01,10)
  sa.analyzeTrans("test2.png",1,0,["vm(2)","vm(3)","vm(4)"],"PULSE(0,1,0.1,0,0,0.1,100.)",0.01,0,3.)
  sa.analyzeManyTrans("test3.png",1,0,"vm(2)",
        [
            "PULSE(0,1,0.1,0,0,0.1,100.)",
            "PULSE(0,1,0.1,0,0,0.3,100.)",
            "PULSE(0,1,0.1,0,0,1.0,100.)",
            "PULSE(0,1,0.1,0,0,3.0,100.)",
        ]
        ,0.01,0,5.)

  with open("example.cir.template") as infile:
    sa = SpiceAnalyzer(infile)
    sa.analyzeAC("test4.png",1,0,["vm(2)","vm(3)","vm(4)"],100,.01,10)
