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
.subckt lowpass 1 3 0
r1 1 2 1
x1 0 2 3 0 opamp
c1 2 3 3e-5
r2 2 3 1
.ends
*
.subckt highpass 1 4 0
r1 1 2 1
c1 2 3 1e-4
x1 0 3 4 0 opamp
r2 3 4 1
.ends
*
.subckt lowpass4 1 4 0
x1 1 2 0 lowpass
x2 2 3 0 lowpass
x3 3 4 0 lowpass
.ends
*
.subckt crrcfltr 1 3 0
x1 1 2 0 highpass
x2 2 3 0 lowpass
.ends
*
.subckt semigaussian 1 3 0
x1 1 2 0 highpass
x2 2 3 0 lowpass4
.ends
*
"""


charge_amp = """charge_amp
*""" + library + """*
*
ci 1 0 100p
rf 1 2 1Meg
cf 1 2 1p
xin 0 1 2 0 opamp
co 2 3 1n
ro 2 3 1K
.end
"""

runs = [
(charge_amp,"Charge Amplifier",["vm(3)"]),
]

for circuit, savename, probes in runs:
  with tempfile.TemporaryFile(mode="w+") as circuitFile:
    circuitFile.write(circuit)
    circuitFile.flush()
    circuitFile.seek(0)
    sa = SpiceAnalyzer(circuitFile)
    sa.analyzeManyTrans(savename+"_trans.png",1,0,probes[0],
          [
              "PULSE(0,1p,10n,0,0,20n,10000u)",
              "PULSE(0,2p,10n,0,0,20n,10000u)",
              "PULSE(0,3p,10n,0,0,20n,10000u)",
          ],
          "1000p",0,"10u",current=True,debug=False)
