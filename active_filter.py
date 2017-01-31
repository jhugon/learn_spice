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
r1 1 2 10K
x1 0 2 3 0 opamp
c1 2 3 10n
r2 2 3 10K
.ends
*
.subckt highpass 1 4 0
r1 1 2 10K
c1 2 3 10n
x1 0 3 4 0 opamp
r2 3 4 10K
.ends
*
.subckt sklowpass2 1 4 0
r1 1 2 10K
r2 2 3 10K
c1 2 4 10n
c2 3 0 10n
x1 3 4 4 0 opamp
.ends
*
.subckt lowpass4 1 4 0
x1 1 2 0 lowpass
x2 2 3 0 lowpass
x3 3 4 0 lowpass
*x4 4 5 0 lowpass
*x5 5 6 0 lowpass
*x6 6 7 0 lowpass
*x7 7 8 0 lowpass
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

active_sk_low_pass2 = """active_sk_low_pass2
*""" + library + """*
*
x1 1 2 0 sklowpass2
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
(active_no_filter,"Active_No_Filter",["3"]),
(active_high_pass,"Active_High_Pass",["2"]),
(active_low_pass,"Active_Low_Pass",["2"]),
(active_sk_low_pass2,"Active_Sallen_Key_Low_Pass_2",["2"]),
(active_crrc_filter,"Active_CRRC_Filter",["2"]),
(active_semigaussian,"Active_Semigaussian_Filter",["2"]),
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
              #"PULSE(0,1,100u,0,0,2000u,10000u)",
              #"PULSE(0,1,100u,0,0,100u,10000u)",
              #"PULSE(0,1,100u,0,0,10u,10000u)",
              "PULSE(0,1,100u,0,0,1u,10000u)",
              "PULSE(0,1,100u,0,0,100n,10000u)",
              "PULSE(0,1,100u,0,0,10n,10000u)",
          ],
          "1u",0,"1000u",debug=False)
