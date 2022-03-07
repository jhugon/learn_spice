#!/usr/bin/env python3

import tempfile
from spiceAnalysis import SpiceAnalyzer
from spiceAnalysis import library as LIBRARY

circuit_template = ""
with open("high_pass.cir.template") as infile:
    circuit_template = infile.read()

sources = [
    "PULSE(0,1,20n,0,0,1p,10000u)",
    "EXP(0,1,20n,1n,20n,5n)",
    "EXP(0,1,20n,1n,20n,25n)",
    "EXP(0,1,20n,1n,20n,50n)",
]

runs = [
    (circuit_template,"differentiator",["2"]),
]

for circuit, savename, probes in runs:
  with tempfile.TemporaryFile(mode="w+") as circuitFile:
    circuitFile.write(circuit)
    circuitFile.flush()
    circuitFile.seek(0)
    sa = SpiceAnalyzer(circuitFile)
    sa.analyzeAC(savename+".png",1,0,probes,1,"1k","1g",debug=False)
    sa.analyzeManyTrans(savename+"_trans_in.png",1,0,"1",
        sources,"500p",0,"200n",debug=False)
    sa.analyzeManyTrans(savename+"_trans.png",1,0,probes[0],
        sources,"500p",0,"200n",debug=False)
