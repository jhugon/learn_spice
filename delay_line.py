#!/usr/bin/env python3

import sys
import tempfile
from spiceAnalysis import SpiceAnalyzer
from spiceAnalysis import library as LIBRARY

circuit_template = ""
with open("../delay_line_lumped_element.cir") as infile:
    for line in infile:
        if line[0] != "J":
            circuit_template += line
circuit_template = circuit_template.replace("Net-_C","").replace("-Pad1_","").replace("/sig_in","98").replace("/sig_out","99").replace("GND","0")
print(circuit_template)
#sys.exit(0)

sources = [
    "PULSE(0,1,50n,0,0,1p,10000u)",
    "EXP(0,1,50n,1n,50n,20n)",
]

runs = [
    (circuit_template,"delay_line",["99"]),
]

for circuit, savename, probes in runs:
  with tempfile.TemporaryFile(mode="w+") as circuitFile:
    circuitFile.write(circuit)
    circuitFile.flush()
    circuitFile.seek(0)
    sa = SpiceAnalyzer(circuitFile)
    sa.analyzeAC(savename+".png","98","0",probes,1,"1k","1g",debug=True)
    sa.analyzeManyTrans(savename+"_trans_in.png",98,0,98,
        sources,"500p",0,"200n",debug=False)
    sa.analyzeManyTrans(savename+"_trans.png",98,0,probes[0],
        sources,"500p",0,"200n",debug=False)
