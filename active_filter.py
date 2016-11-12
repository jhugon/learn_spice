#!/usr/bin/env python3

import tempfile
from spiceAnalysis import SpiceAnalyzer

active_no_filter = """active_no_filter
*
* 1 input+, 2 input-, 3 output, 4 ground
.subckt opamp 1 2 3 4
E1 3 4 1 2 1e8
.ends
*
r1 1 2 1K
x1 0 2 3 0 opamp
r2 2 3 1K
.end
"""


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
c1 2 3 10n
r2 2 3 1K
.end
"""

active_band_pass = """active_band_pass
*
* 1 input+, 2 input-, 3 output, 4 ground
.subckt opamp 1 2 3 4
E1 3 4 1 2 1e8
.ends
*
* first high pass
r1 1 2 1K
c1 2 3 10n
x1 0 3 4 0 opamp
r2 3 4 2K
* now low pass
r3 4 5 1K
x2 0 5 6 0 opamp
c2 5 6 10n
r4 5 6 1K
.end
"""

active_band_pass2 = """active_band_pass2
*
* 1 input+, 2 input-, 3 output, 4 ground
.subckt opamp 1 2 3 4
E1 3 4 1 2 1e8
.ends
*
* first low pass
r1 1 2 1K
x1 0 2 3 0 opamp
c1 2 3 10n
r2 2 3 1K
*
r3 3 4 1K
c2 4 5 100n
x2 0 5 6 0 opamp
r4 5 6 2K
.end
"""

runs = [
(active_no_filter,"Active_No_Filter",["vm(3)"]),
(active_high_pass,"Active_High_Pass",["vm(4)"]),
(active_low_pass,"Active_Low_Pass",["vm(3)"]),
(active_band_pass,"Active_band_Pass",["vm(6)"]),
(active_band_pass2,"Active_band_Pass2",["vm(6)"]),
]

for circuit, savename, probes in runs:
  with tempfile.TemporaryFile(mode="w+") as circuitFile:
    circuitFile.write(circuit)
    circuitFile.flush()
    circuitFile.seek(0)
    sa = SpiceAnalyzer(circuitFile)
    sa.analyzeAC(savename+".png",1,0,probes,1,"1k","1000k",debug=False)
    sa.analyzeManyTrans(savename+"_trans.png",1,0,"vm(3)",
          [
              "PULSE(0,1,10u,1u,0,2000u,10000u)",
              "PULSE(0,1,10u,100u,0,2000u,10000u)",
              "PULSE(0,1,10u,500u,0,2000u,10000u)",
              "PULSE(0,1,500u,1u,0,200u,10000u)",
              "PULSE(0,1,500u,0,0,10u,10000u)",
          ],
          "1u",0,"1000u",debug=False)
