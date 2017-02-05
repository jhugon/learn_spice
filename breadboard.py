#!/usr/bin/env python3

import sys
import tempfile
import re
from spiceAnalysis import SpiceAnalyzer

library = ""

with open("LM324.5_1") as lm324file:
  lm324str = lm324file.read()
  library += lm324str[:-1]
  # 1.2MHz bandwidth quad op-amp single/dual supply PDIP
  #* (REV N/A)      SUPPLY VOLTAGE: 5V
  #* CONNECTIONS:   NON-INVERTING INPUT
  #*                | INVERTING INPUT
  #*                | | POSITIVE POWER SUPPLY
  #*                | | | NEGATIVE POWER SUPPLY
  #*                | | | | OUTPUT
  #*                | | | | |
  #.SUBCKT LM324    1 2 3 4 5

with open("LM318.301") as chipfile:
  chipstr = chipfile.read()
  library += chipstr[:-1]
  # 15 MHz bandwidth high slew rate Double PDIP
  #* (REV N/A)      SUPPLY VOLTAGE: +/-15V
  #* CONNECTIONS:   NON-INVERTING INPUT
  #*                | INVERTING INPUT
  #*                | | POSITIVE POWER SUPPLY
  #*                | | | NEGATIVE POWER SUPPLY
  #*                | | | | OUTPUT
  #*                | | | | |
  #.SUBCKT LM318    1 2 3 4 5

with open("AD8651.cir") as chipfile:
  chipstr = chipfile.read()
  library += chipstr
  # 50 MHz precision low-noise amp, smd
  #*			noninverting input
  #*			|	inverting input
  #*			|	|	 positive supply
  #*			|	|	 |	 negative supply
  #*			|	|	 |	 |	 output
  #*			|	|	 |	 |	 |
  #*			|	|	 |	 |	 |
  #.SUBCKT AD8651		1	2	99	50	45

with open("ada4627.cir") as chipfile:
  chipstr = chipfile.read()
  library += chipstr
  # 19 MHz SMD JFET input 
  #*                 non-inverting input
  #*                 |   inverting input
  #*                 |   |    positive supply
  #*                 |   |    |    negative supply
  #*                 |   |    |    |    output
  #*                 |   |    |    |    |
  #.SUBCKT ADA4627   1   2   99   50   30

with open("ada4637.cir") as chipfile:
  chipstr = chipfile.read()
  library += chipstr
  # 79 MHz SMD JFET input 
  #*            	  non-inverting input
  #*                 |   inverting input
  #*                 |   |    positive supply
  #*                 |   |    |    negative supply
  #*                 |   |    |    |    output
  #*                 |   |    |    |    |
  #.SUBCKT ADA4637   1   2   99   50   30

##enable this if you want to use the JFET chips
#library += """*
#    .OPTIONS GMIN=0.01p
#    .OPTIONS ABSTOL=0.01pA
#    .OPTIONS ITL1=500
#    .OPTIONS ITL2=200
#    .OPTIONS ITL4=100
#"""

with open("TLC2274.101") as chipfile:
  chipstr = chipfile.read()
  library += chipstr[:-1]
  #CMOS Quad RR op amp DIP 2 MHz
  #* CONNECTIONS:   NON-INVERTING INPUT
  #*                | INVERTING INPUT
  #*                | | POSITIVE POWER SUPPLY
  #*                | | | NEGATIVE POWER SUPPLY
  #*                | | | | OUTPUT
  #*                | | | | |
  #.SUBCKT TLC2274  1 2 3 4 5

with open("OPA277.txt") as chipfile:
  chipstr = chipfile.read()
  newchipstr = """* comatibility stuff
.func LIMIT(x,a,b) {min(max(x, a), b)}
.func PWR(x,a) {abs(x) ** a}
.func PWRS(x,a) {sgn(x) * PWR(x,a)}
.func stp(x) {u(x)}
  """
  for line in chipstr.split("\n"):
    vswitchmodelMatch = re.search("^(\.model.*)vswitch",line.lower())
    if vswitchmodelMatch:
      line = line.lower()
      ronMatch = re.search(r"\.model.+ron=([a-zA-Z0-9_-]+)",line)
      if not ronMatch:
        print("In VSWITCH Converter, couldn't find RON. Exiting.")
        sys.exit(1)
      roffMatch = re.search(r"\.model.+roff=([a-zA-Z0-9_-]+)",line)
      if not ronMatch:
        print("In VSWITCH Converter, couldn't find ROFF. Exiting.")
        sys.exit(1)
      vonMatch = re.search(r"\.model.+von=([a-zA-Z0-9_-]+)",line)
      if not ronMatch:
        print("In VSWITCH Converter, couldn't find VON. Exiting.")
        sys.exit(1)
      voffMatch = re.search(r"\.model.+voff=([a-zA-Z0-9_-]+)",line)
      if not ronMatch:
        print("In VSWITCH Converter, couldn't find VOFF. Exiting.")
        sys.exit(1)
      voff = voffMatch.group(1)
      von = vonMatch.group(1)
      roff = roffMatch.group(1)
      ron = ronMatch.group(1)
      newline = "{}aswitch(cntl_off={} cntl_on={} r_off={} r_on={} log=TRUE)".format(vswitchmodelMatch.group(1),voff,von,roff,ron)
      newchipstr += newline + "\n"
    else:
      newchipstr += line + "\n"
  library += newchipstr
  # 1 MHz DIP High Precision
  #.SUBCKT OPA277 +IN -IN V+ V- Vout

with open("OPA227.MOD") as chipfile:
  chipstr = chipfile.read()
  library += chipstr
  # 8 MHz DIP Precision Low Noise
  #* PINOUT        3   2   7  4  6
  #* PINOUT ORDER +IN -IN +V -V OUT
  #.SUBCKT OPA227 3 2 7 4 6


library += """

*
* 1 input+, 2 input-, 3 output, 4 +V supply, 5 -V supply, 6 ground
.subckt opamp 1 2 3 4 5 6
*E1 3 6 1 2 1e8
****x1 1 2 4 5 3 LM318
****x1 1 2 4 5 3 AD8651
****x1 1 2 4 5 3 ADA4627
****x1 1 2 4 5 3 ADA4637
*
x1 1 2 4 5 3 LM324
*x1 1 2 4 5 3 TLC2274
*x1 1 2 4 5 3 OPA277
*x1 1 2 4 5 3 OPA227
.ends
*
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
.subckt lowpass4 1 5 10 11 0
x1 1 2 10 11 0 lowpass
x2 2 3 10 11 0 lowpass
x3 3 4 10 11 0 lowpass
x4 4 5 10 11 0 lowpass
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

active_low_pass4 = """active_low_pass4
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
x1 1 2 99 100 0 lowpass4
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
(active_low_pass4,"Bread_Low_Pass4",["2"]),
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
