#!/usr/bin/env python3

import re
import os
import tempfile
import subprocess
import numpy
import matplotlib.pyplot as mpl

def getData(logfile):
  xtitle = None
  ytitle = None
  xdata = []
  ydata = []
  foundStart1 = False
  foundStart2 = False
  foundStart3 = False
  for line in logfile:
    if foundStart3:
      datamatch = re.match(r"^(\d+)\s+([0-9eE.+-]+)\s+([0-9eE.+-]+)\s*$",line)
      if datamatch:
        x = datamatch.group(2)
        y = datamatch.group(3)
        x = float(x)
        y = float(y)
        xdata.append(x)
        ydata.append(y)
    elif foundStart2:
      if re.match(r"^-+\s*$",line):
          foundStart3 = True
    elif foundStart1:
      start2match = re.match(r"^Index\s+(\w+)\s+([a-zA-Z0-9\(\)]+)\s*$",line)
      if start2match:
        xtitle = start2match.group(1)
        ytitle = start2match.group(2)
        foundStart2 = True
    elif re.match(r"^-+\s*$",line):
      foundStart1 = True
  xdata = numpy.array(xdata)
  ydata = numpy.array(ydata)
  return xdata, ydata, xtitle, ytitle

def plotData(filename):
  xdata,ydata, xtitle,ytitle = getData(filename)
  fig, ax = mpl.subplots()
  ax.plot(xdata,ydata)
  ax.set_xlabel(xtitle)
  ax.set_ylabel(ytitle)
  fig.savefig(filename+".png")

class TemplateModifier(object):

  def __init__(self,circuitTemplateFilename):
    self.templateFilename = circuitTemplateFilename
    self.sourceLines = []
    self.analyzerLines = [".options noacct"]
    self.tmpFiles = []

  def addACSource(self,name,nodePlus,nodeMinus,mag,phase=0.):
    line = "{0} {1:d} {2:d} AC {3} {4}".format(name,nodePlus,nodeMinus,mag,phase)
    self.sourceLines.append(line)

  def addACAnalysis(self,nodePlus,nodeMinus,pointsPerDecade,fstart,fstop):
    line = ".ac dec {0:d} {1} {2}".format(pointsPerDecade,fstart,fstop)
    self.sourceLines.append(line)
    line = ".print ac vm({0:d},{1:d})".format(nodePlus,nodeMinus)
    self.sourceLines.append(line)

  def getFile(self):
    result = tempfile.NamedTemporaryFile(mode="w+")
    with open(self.templateFilename) as infile:
        for line in infile:
          if line[:4] != ".end":
            result.write(line)
    for line in self.sourceLines:
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

def analyzeAC(circuitTemplate,outFileName,inNodePlus,inNodeMinus,outNodePlus,outNodeMinus,mag,fstart,fstop,pointsPerDecade=5):
  template = TemplateModifier(circuitTemplate)
  template.addACSource("vac",inNodePlus,inNodeMinus,mag)
  template.addACAnalysis(outNodePlus,outNodeMinus,pointsPerDecade,fstart,fstop)
  circuitFileName = template.getFile()
  if True:
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
    
  if True:
    output.seek(0)
    print(output.read())
  output.seek(0)
  xdata,ydata, xtitle, ytitle = getData(output)
  ydataDB = 20*numpy.log(ydata/mag)
  fig, ax = mpl.subplots()
  ax.plot(xdata,ydataDB)
  ax.set_xlabel(xtitle)
  ax.set_ylabel("V$_{{{},{}}}$ [dBc]".format(outNodePlus,outNodeMinus))
  fig.savefig(outFileName)
  
analyzeAC("example.cir.template","test.png",1,0,2,0,100,.01,10)
