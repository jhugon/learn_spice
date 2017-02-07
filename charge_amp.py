#!/usr/bin/env python3

import tempfile
from spiceAnalysis import SpiceAnalyzer

library = """*
* 1 input+, 2 input-, 3 output, 4 ground
.subckt opamp 1 2 3 4
E1 3 4 1 2 1e8
.ends
*
* all filters below are input, output, grnd
*
*
"""

charge_amp = """charge_amp
*""" + library + """*
*
ri 1 0 100
ci 1 0 1p
rf 1 2 1Meg
cf 1 2 1n
xin 0 1 2 0 opamp
rc 2 3 1Meg
cc 2 3 1n
rt 3 4 100k
ct 3 4 1n
xt 0 3 4 0 opamp
.end
"""

runs = [
(charge_amp,"Charge Amplifier",["3"]),
]

for circuit, savename, probes in runs:
  with tempfile.TemporaryFile(mode="w+") as circuitFile:
    circuitFile.write(circuit)
    circuitFile.flush()
    circuitFile.seek(0)
    sa = SpiceAnalyzer(circuitFile)
    #sa.analyzeAC(savename+"_ac.png",1,0,probes,1,"10","100M",debug=False)
    sa.analyzeManyTrans(savename+"_singletrans.png",1,0,probes[0],
          [
              "PULSE(0,1n,0,0,0,10u,1)",
              "PULSE(0,2n,0,0,0,10u,1)",
              "PULSE(0,3n,0,0,0,10u,1)",
          ],
          "100n",0,"1m",current=True,debug=False)
    #sa.analyzeManyTrans(savename+"_manytrans.png",1,0,probes[0],
    #      [
    #          "PULSE(0,1n,0,0,0,10u,10u)",
    #          "PULSE(0,1n,0,0,0,10u,20u)",
    #          "PULSE(0,1n,0,0,0,10u,30u)",
    #          "PULSE(0,1n,0,0,0,10u,40u)",
    #      ],
    #      "100n",0,"10m",current=True,debug=False)
