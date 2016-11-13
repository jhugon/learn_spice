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

active_no_filter = """active_no_filter
*""" + library + """*
*
r1 1 2 1K
x1 0 2 3 0 opamp
r2 2 3 1K
.end
"""


active_high_pass = """active_high_pass
*""" + library + """*
*
x1 1 2 0 highpass
.end
"""

active_low_pass = """active_low_pass
*""" + library + """*
*
x1 1 2 0 lowpass4
.end
"""

active_crrc_filter = """active_crrc_filter
*""" + library + """*
*
x1 1 2 0 crrcfltr
.end
"""

active_semigaussian = """active_semigaussian
*""" + library + """*
*
x1 1 2 0 semigaussian
.end
"""

runs = [
(active_no_filter,"Active_No_Filter",["vm(3)"]),
(active_high_pass,"Active_High_Pass",["vm(2)"]),
(active_low_pass,"Active_Low_Pass",["vm(2)"]),
(active_crrc_filter,"Active_CRRC_Filter",["vm(2)"]),
(active_semigaussian,"Active_Semigaussian_Filter",["vm(2)"]),
]

for circuit, savename, probes in runs:
  with tempfile.TemporaryFile(mode="w+") as circuitFile:
    circuitFile.write(circuit)
    circuitFile.flush()
    circuitFile.seek(0)
    sa = SpiceAnalyzer(circuitFile)
    sa.analyzeAC(savename+".png",1,0,probes,1,"1","10000k",debug=False)
    sa.analyzeManyTrans(savename+"_trans.png",1,0,probes[0],
          [
              "PULSE(0,1,100u,0,0,2000u,10000u)",
              "PULSE(0,1,100u,0,0,100u,10000u)",
          ],
          "1u",0,"1000u",debug=False)