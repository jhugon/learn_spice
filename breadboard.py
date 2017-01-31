#!/usr/bin/env python3

import sys
import tempfile
from spiceAnalysis import SpiceAnalyzer

library = ""

with open("LM324.5_1") as lm324file:
  lm324str = lm324file.read()
  library += lm324str[:-1]
  #* (REV N/A)      SUPPLY VOLTAGE: 5V
  #* CONNECTIONS:   NON-INVERTING INPUT
  #*                | INVERTING INPUT
  #*                | | POSITIVE POWER SUPPLY
  #*                | | | NEGATIVE POWER SUPPLY
  #*                | | | | OUTPUT
  #*                | | | | |
  #.SUBCKT LM324    1 2 3 4 5

library += """

*
* 1 input+, 2 input-, 3 output, 4 ground
* 1 input+, 2 input-, 3 +V supply, 4 -V supply, 5 output
.subckt opamp 1 2 3 4 5
E1 5 4 1 2 1e8
.ends
*
* all filters below are input, output, grnd
* all filters below are input, output, grnd, +supply
*
.subckt lowpass 1 3 0 10
r1 1 2 10K
x1 0 2 10 0 3 LM324
c1 2 3 10n
r2 2 3 10K
.ends
*
.subckt highpass 1 4 0 10
r1 1 2 10K
c1 2 3 10n
x1 0 3 10 0 4 LM324
r2 3 4 10K
.ends
*
.subckt sklowpass2 1 4 0 10
r1 1 2 10K
r2 2 3 10K
c1 2 4 10n
c2 3 0 10n
x1 3 4 10 0 4 LM324
.ends
*
.subckt lowpass4 1 4 0 10
x1 1 2 0 10 lowpass
x2 2 3 0 10 lowpass
x3 3 4 0 10 lowpass
*x4 4 5 0 10 lowpass
*x5 5 6 0 10 lowpass
*x6 6 7 0 10 lowpass
*x7 7 8 0 10 lowpass
.ends
*
.subckt crrcfltr 1 3 0 10
x1 1 2 0 10 highpass
x2 2 3 0 10 lowpass
.ends
*
.subckt semigaussian 1 3 0 10
x1 1 2 0 10 highpass
x2 2 3 0 10 lowpass4
.ends
*
"""

active_high_pass = """active_high_pass
*""" + library + """*
*
v99 99 0 DC 5
x1 1 2 0 99 highpass
.end
"""

active_low_pass = """active_low_pass
*""" + library + """*
*
v99 99 0 DC 5
x1 1 2 0 99 lowpass4
.end
"""

active_sk_low_pass2 = """active_sk_low_pass2
*""" + library + """*
*
v99 99 0 DC 5
x1 1 2 0 99 sklowpass2
.end
"""

active_crrc_filter = """active_crrc_filter
*""" + library + """*
*
v99 99 0 DC 5
x1 1 2 0 99 crrcfltr
.end
"""

active_semigaussian = """active_semigaussian
*""" + library + """*
*
v99 99 0 DC 5
x1 1 2 0 99 semigaussian
.end
"""

runs = [
#(active_no_filter,"Bread_No_Filter",["3"]),
(active_high_pass,"Bread_High_Pass",["2"]),
(active_low_pass,"Bread_Low_Pass",["2"]),
(active_sk_low_pass2,"Bread_Sallen_Key_Low_Pass_2",["2"]),
(active_crrc_filter,"Bread_CRRC_Filter",["2"]),
(active_semigaussian,"Bread_Semigaussian_Filter",["2"]),
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
              #"PULSE(0,1,100u,0,0,1u,10000u)",
              #"PULSE(0,1,100u,0,0,100n,10000u)",
              #"PULSE(0,1,100u,0,0,10n,10000u)",
               "SIN(0,1,1000,0,0)", # offset, amp, freq, delay, damping
          ],
          "10u",0,"10000u",debug=False)
