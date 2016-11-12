#!/usr/bin/env python3

import tempfile
from spiceAnalysis import SpiceAnalyzer

active_high_pass = """active_high_pass
*
* 1 input+, 2 input-, 3 output, 4 ground
.subckt opamp 1 2 3 4
E1 3 4 1 2 1e8
.ends
*
r1 1 2 1K
c1 2 3 10n
x1 0 3 4 0 opamp
r2 3 4 1K
.end
"""

active_low_pass = """active_low_pass
*
* 1 input+, 2 input-, 3 output, 4 ground
.subckt opamp 1 2 3 4
E1 3 4 1 2 1e8
.ends
*
r1 1 2 1K
x1 0 2 3 0 opamp
c1 2 3 1n
r2 2 3 1K
.end
"""

for circuit, savename, probes in [(active_high_pass,"Active_High_Pass",["vm(4)"]),(active_low_pass,"Active_Low_Pass",["vm(3)"])]:
  with tempfile.TemporaryFile(mode="w+") as circuitFile:
    circuitFile.write(circuit)
    circuitFile.flush()
    circuitFile.seek(0)
    sa = SpiceAnalyzer(circuitFile)
    sa.analyzeAC(savename+".png",1,0,probes,1,"1k","1000k",debug=False)
    #sa.analyzeManyTrans("HighLow2.png",1,0,"vm(3)",
    #      [
    #          "PULSE(0,1,10u,0,0,1u,100u)",
    #          "PULSE(0,1,10u,0,0,10u,100u)",
    #          "PULSE(0,1,10u,0,0,30u,100u)",
    #          "PULSE(0,1,10u,0,0,200u,100u)",
    #      ],
    #      "0.25u",0,"100u",debug=False)
