#!/usr/bin/env python3

import sys
import tempfile
import re
from spiceAnalysis import SpiceAnalyzer
from spiceAnalysis import library as LIBRARY

library = LIBRARY

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

#with open("TLC2274.101") as chipfile:
#  chipstr = chipfile.read()
#  library += chipstr[:-1]
#  #CMOS Quad RR op amp DIP 2 MHz
#  #* CONNECTIONS:   NON-INVERTING INPUT
#  #*                | INVERTING INPUT
#  #*                | | POSITIVE POWER SUPPLY
#  #*                | | | NEGATIVE POWER SUPPLY
#  #*                | | | | OUTPUT
#  #*                | | | | |
#  #.SUBCKT TLC2274  1 2 3 4 5

#with open("OPA277.txt") as chipfile:
#  chipstr = chipfile.read()
#  newchipstr = """* comatibility stuff
#.func LIMIT(x,a,b) {min(max(x, a), b)}
#.func PWR(x,a) {abs(x) ** a}
#.func PWRS(x,a) {sgn(x) * PWR(x,a)}
#.func stp(x) {u(x)}
#  """
#  for line in chipstr.split("\n"):
#    vswitchmodelMatch = re.search("^(\.model.*)vswitch",line.lower())
#    if vswitchmodelMatch:
#      line = line.lower()
#      ronMatch = re.search(r"\.model.+ron=([a-zA-Z0-9_-]+)",line)
#      if not ronMatch:
#        print("In VSWITCH Converter, couldn't find RON. Exiting.")
#        sys.exit(1)
#      roffMatch = re.search(r"\.model.+roff=([a-zA-Z0-9_-]+)",line)
#      if not ronMatch:
#        print("In VSWITCH Converter, couldn't find ROFF. Exiting.")
#        sys.exit(1)
#      vonMatch = re.search(r"\.model.+von=([a-zA-Z0-9_-]+)",line)
#      if not ronMatch:
#        print("In VSWITCH Converter, couldn't find VON. Exiting.")
#        sys.exit(1)
#      voffMatch = re.search(r"\.model.+voff=([a-zA-Z0-9_-]+)",line)
#      if not ronMatch:
#        print("In VSWITCH Converter, couldn't find VOFF. Exiting.")
#        sys.exit(1)
#      voff = voffMatch.group(1)
#      von = vonMatch.group(1)
#      roff = roffMatch.group(1)
#      ron = ronMatch.group(1)
#      newline = "{}aswitch(cntl_off={} cntl_on={} r_off={} r_on={} log=TRUE)".format(vswitchmodelMatch.group(1),voff,von,roff,ron)
#      newchipstr += newline + "\n"
#    else:
#      newchipstr += line + "\n"
#  library += newchipstr
#  # 1 MHz DIP High Precision
#  #.SUBCKT OPA277 +IN -IN V+ V- Vout
#
#with open("OPA227.MOD") as chipfile:
#  chipstr = chipfile.read()
#  library += chipstr
#  # 8 MHz DIP Precision Low Noise
#  #* PINOUT        3   2   7  4  6
#  #* PINOUT ORDER +IN -IN +V -V OUT
#  #.SUBCKT OPA227 3 2 7 4 6


with open("ada4001.cir") as chipfile:
  chipstr = chipfile.read()
  library += chipstr
  # 16.7 MHz SMD JFET, Low Noise, Low Ib, RRO
  #*                non-inverting input
  #*                | inverting input
  #*                | | positive supply
  #*                | | |  negative supply
  #*                | | |  |  output
  #*                | | |  |  |
  #.SUBCKT ADA4001  1 2 99 50 45


library += """

*
* 1 input+, 2 input-, 3 output, 4 +V supply, 5 -V supply, 6 ground
.subckt opamp 1 2 3 4 5 6
****x1 1 2 4 5 3 LM318
****x2 1 2 4 5 3 AD8651
****x3 1 2 4 5 3 ADA4627
****x4 1 2 4 5 3 ADA4637
*
*x0 1 2 3 4 5 6 idealopamp
*x5 1 2 4 5 3 LM324
*x6 1 2 4 5 3 TLC2274
*x7 1 2 4 5 3 OPA277
*x8 1 2 4 5 3 OPA227
*x9 1 2 3 4 5 6 idealopamp1k
x11 1 2 4 5 3 ADA4001
.ends
*
*
.subckt buffer 1 0 2 99 100 0
x1 1 2 2 99 100 0 opamp
.ends
*
.subckt semigausmfb50 1 0 3 99 100 0
r1 1 2 6.771912e2
c1 2 0 1n
x1 2 0 3 99 100 0 buffer
.ends
*
.subckt semigausmfb51 1 0 4 99 100 0
r2 1 2 4.8093e2
r1 2 3 34.135
r3 2 4 4.8093e2
c1 2 0 2.4n
c2 3 4 1n
x1 0 3 4 99 100 0 opamp
.ends
*
.subckt semigausmfb52 1 0 4 99 100 0
r2 1 2 2.9288e2
r1 2 3 23.69
r3 2 4 2.9288e2
c1 2 0 3.9n
c2 3 4 1n
x1 0 3 4 99 100 0 opamp
.ends
*
*
"""

active_mfb = """active_mfb
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
x2 1 0 2 99 100 0 semigausmfb50
.end
"""

active_semigaussian = """active_semigaussian
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
*c1 1 3 10n
*r1 3 0 10k
*xb1 3 0 4 99 100 0 buffer
xb1 1 0 4 99 100 0 buffer
x0 4 0 5 99 100 0 semigausmfb50
x1 5 0 6 99 100 0 semigausmfb51
x2 6 0 7 99 100 0 semigausmfb52
.end
"""

active_crrc4 = """active_crrc4
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
c1 1 3 1n
r1 3 0 1.58k
xb1 3 0 4 99 100 0 buffer
r2 3 4 1.58k
c2 4 0 1n
xb2 4 0 5 99 100 0 buffer
r3 5 6 1.58k
c3 6 0 1n
xb3 6 0 7 99 100 0 buffer
r4 7 8 1.58k
c4 8 0 1n
xb4 8 0 9 99 100 0 buffer
r5 9 10 1.58k
c5 10 0 1n
xb5 10 0 11 99 100 0 buffer
.end
"""

runs = [
(active_mfb,"MFB_Filter",["2"]),
(active_semigaussian,"Semigaussian_Filter",["7"]),
(active_crrc4,"CRRC4_Filter",["11"]),
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
              "PULSE(0,1,0,0,0,1u,1)",
              "PULSE(0,1,0,0,0,2.5u,1)",
              "PULSE(0,1,0,0,0,5u,1)",
              "PULSE(0,1,0,0,0,7.5u,1)",
              "PULSE(0,1,0,0,0,10u,1)",
              "PULSE(0,1,0,0,0,100u,1)",
              # "SIN(0,1,500,0,0)", # offset, amp, freq, delay, damping
          ],
          "0.1u",0,"40u",debug=False)
    sa.analyzeManyTrans(savename+"_trans2.png",1,0,probes[0],
          [
            
              # pulse args are initial val, pulsed val, delay, rise time, fall time, pulse width, period.
              "PULSE(0,1,0,0,0,2.5u,7.5u)",
              "PULSE(0,1,0,0,0,2.5u,5.0u)",
              "PULSE(0,1,0,0,0,2.5u,4.0u)",
              "PULSE(0,1,0,0,0,2.5u,3.0u)",
              # "SIN(0,1,500,0,0)", # offset, amp, freq, delay, damping
          ],
          "0.1u",0,"50u",debug=False)
