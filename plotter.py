#!/usr/bin/env python3

import re
import os
import tempfile
import subprocess
import numpy
import matplotlib.pyplot as mpl

def getData(logfile,nYCols):
  xtitle = None
  ytitles = []
  xdata = []
  ydata = []
  for iCol in range(nYCols):
    ydata.append([])
  foundStart1 = False
  foundStart2 = False
  foundStart3 = False
  for line in logfile:
    if foundStart3:
      reStr = r"^(\d+)\s+([0-9eE.+-]+)"+("\s+([0-9eE.+-]+)")*nYCols+"\s*$"
      datamatch = re.match(reStr,line)
      if datamatch:
        x = datamatch.group(2)
        x = float(x)
        xdata.append(x)
        for iCol in range(nYCols):
          y = datamatch.group(3+iCol)
          y = float(y)
          ydata[iCol].append(y)
    elif foundStart2:
      if re.match(r"^-+\s*$",line):
          foundStart3 = True
    elif foundStart1:
      reStr = r"^Index\s+(\w+)"+nYCols*("\s+([a-zA-Z0-9\(\)]+)")+"\s*$"
      start2match = re.match(reStr,line)
      if start2match:
        xtitle = start2match.group(1)
        for iCol in range(nYCols):
          ytitles.append(start2match.group(2+iCol))
        foundStart2 = True
    elif re.match(r"^-+\s*$",line):
      foundStart1 = True
  xdata = numpy.array(xdata)
  ydata = numpy.array(ydata)
  return xdata, ydata, xtitle, ytitles

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

  def addACAnalysis(self,nodesToPrint,pointsPerDecade,fstart,fstop):
    line = ".ac dec {0:d} {1} {2}".format(pointsPerDecade,fstart,fstop)
    self.sourceLines.append(line)
    line = ".print ac"
    for node in nodesToPrint:
      if type(node) == list:
        line += " vm({0:d},{1:d})".format(node[0],node[1])
      else:
        line += " vm({0:d})".format(node)
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

def analyzeAC(circuitTemplate,outFileName,inNodePlus,inNodeMinus,outNodes,mag,fstart,fstop,pointsPerDecade=5):
  template = TemplateModifier(circuitTemplate)
  template.addACSource("vac",inNodePlus,inNodeMinus,mag)
  template.addACAnalysis(outNodes,pointsPerDecade,fstart,fstop)
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
  xdata,ydata, xtitle, ytitle = getData(output,len(outNodes))
  if True:
    print(xdata)
    for iCol in range(len(outNodes)):
        print(ydata[iCol])
  ydataDB = 20*numpy.log(ydata/mag)
  fig, ax = mpl.subplots()
  for iCol in range(len(outNodes)):
    outNode = outNodes[iCol]
    label = None
    if type(outNode) == list:
      label = "Node {} to {}".format(outNode[0],outNode[1])
    else:
      label = "Node {} to 0".format(outNode)
    ax.semilogx(xdata,ydataDB[iCol],label=label)
  ax.set_xlabel(xtitle)
  ax.set_ylabel("V [dBc]")
  ax.legend()
  fig.savefig(outFileName)
  
analyzeAC("example.cir.template","test.png",1,0,[2,3],100,.01,10)
