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
* 1 input+, 2 input-, 3 +V supply, 4 -V supply, 5 output
* 1 input+, 2 input-, 3 output, 4 +V supply, 5 -V supply, 6 ground
.subckt opamp 1 2 3 4 5 6
*E1 3 6 1 2 1e8
x1 1 2 4 5 3 LM324
.ends
*
* all filters below are input, output, grnd, +supply
* all filters below are input, output, +supply, -supply, ground
*
.subckt lowpass 1 3 10 11 0
r1 1 2 1.58K
x1 0 2 3 10 11 0 opamp
c1 2 3 1n
r2 2 3 1.58K
.ends
*
.subckt highpass 1 4 10 11 0
r1 1 2 1.58K
c1 2 3 1n
x1 0 3 4 10 11 0 opamp
r2 3 4 1.58K
.ends
*
.subckt sklowpass2 1 4 10 11 0
r1 1 2 1.58K
r2 2 3 1.58K
c1 2 4 1n
c2 3 0 1n
x1 3 4 4 10 11 0 opamp
.ends
*
.subckt lowpass4 1 4 10 11 0
x1 1 2 10 11 0 lowpass
x2 2 3 10 11 0 lowpass
x3 3 4 10 11 0 lowpass
*x4 4 5 10 11 0 lowpass
*x5 5 6 10 11 0 lowpass
*x6 6 7 10 11 0 lowpass
*x7 7 8 10 11 0 lowpass
.ends
*
.subckt crrcfltr 1 3 10 11 0
x1 1 2 10 11 0 highpass
x2 2 3 10 11 0 lowpass
.ends
*
.subckt semigaussian 1 3 10 11 0
x1 1 2 10 11 0 highpass
x2 2 3 10 11 0 lowpass4
.ends
*
"""

active_emitter_follower = """active_emitter_follower
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
x1 1 3 3 99 100 0 opamp
.end
"""

active_10gain = """active_10gain
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
r1 1 2 300
r2 2 3 3k
x1 0 2 3 99 100 0 opamp
.end
"""

active_100gain = """active_100gain
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
r1 1 2 100
r2 2 3 10k
x1 0 2 3 99 100 0 opamp
.end
"""

active_high_pass = """active_high_pass
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
x1 1 2 99 100 0 highpass
.end
"""

active_low_pass = """active_low_pass
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
x1 1 2 99 100 0 lowpass
.end
"""

active_sk_low_pass2 = """active_sk_low_pass2
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
x1 1 2 99 100 0 sklowpass2
.end
"""

active_crrc_filter = """active_crrc_filter
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
x1 1 2 99 100 0 crrcfltr
.end
"""

active_semigaussian = """active_semigaussian
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
x1 1 2 99 100 0 semigaussian
.end
"""

runs = [
(active_emitter_follower,"Bread_Emitter_follower",["3"]),
(active_10gain,"Bread_10gain",["3"]),
(active_100gain,"Bread_100gain",["3"]),
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
    sa.analyzeAC(savename+".png",1,0,probes,1,"100","100Meg",debug=False)
    sa.analyzeManyTrans(savename+"_trans.png",1,0,probes[0],
          [
            
              # pulse args are initial val, pulsed val, delay, rise time, fall time, pulse width, period.
              "PULSE(0,1,10u,0,0,2000u,10000u)",
              # "SIN(0,1,500,0,0)", # offset, amp, freq, delay, damping
          ],
          "1u",0,"100u",debug=False)
