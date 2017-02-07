#!/usr/bin/env python3

import tempfile
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

with open("OPA227.MOD") as chipfile:
  chipstr = chipfile.read()
  library += chipstr
  # 8 MHz DIP Precision Low Noise
  #* PINOUT        3   2   7  4  6
  #* PINOUT ORDER +IN -IN +V -V OUT
  #.SUBCKT OPA227 3 2 7 4 6


library += """*
* 1 input+, 2 input-, 3 output, 4 +V supply, 5 -V supply, 6 ground
.subckt opamp 1 2 3 99 100 4
x0 1 2 3 99 100 4 idealopamp
*x1 1 2 99 100 3 LM324
*x2 1 2 99 100 3 LM318
*x3 1 2 99 100 3 ADA4627
.ends
*
*
*
* 1e-4 wT, 1e-3 wF, 1e5 gain
* 1 input+, 2 input-, 3 output, 4 +V supply, 5 -V supply, 6 ground
.subckt charge_amp 2 0 5 99 100 0
rf 2 3 100k
cf 2 3 10n
xin 0 2 3 99 100 0 opamp
rc 3 4 1Meg
cc 3 4 1n
rt 4 5 10k
ct 4 5 10n
xt 0 4 5 99 100 0 opamp
.ends
"""

charge_amp = """charge_amp
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
ri 1 0 1Meg
ci 1 0 100p
* blocking cap or not
cb 1 2 1u
*rb 1 2 0
*
xca 2 0 3 99 100 0 charge_amp
.end
"""

charge_amp_test_cap = """charge_amp_test_cap
*""" + library + """*
*
v99 99 0 DC 5
v100 100 0 DC -5
ct 1 2 100n
rt 1 2 1000G
*
xca 2 0 3 99 100 0 charge_amp
.end
"""

runs = [
(charge_amp,"Charge_Amplifier",["3"]),
]

for circuit, savename, probes in runs:
  with tempfile.TemporaryFile(mode="w+") as circuitFile:
    circuitFile.write(circuit)
    circuitFile.flush()
    circuitFile.seek(0)
    sa = SpiceAnalyzer(circuitFile)
    #sa.analyzeAC(savename+"_ac.png",1,0,probes,1,"10","100M",current=True,debug=False)
    sa.analyzeManyTrans(savename+"_singletrans.png",1,0,probes[0],
          [
              "PULSE(0,1n,100u,0,0,1u,1)",
              "PULSE(0,1n,100u,0,0,2u,1)",
              "PULSE(0,2n,100u,0,0,2u,1)",
          ],
          "100n",0,"1m",current=True,debug=False)
    sa.analyzeManyTrans(savename+"_manytrans.png",1,0,probes[0],
          [
              "PULSE(0,1n,0,0,0,1u,300u)",
              "PULSE(0,1n,0,0,0,1u,500u)",
              "PULSE(0,1n,0,0,0,1u,700u)",
          ],
          "100n",0,"1.5m",current=True,debug=False)

## Test cap time
with tempfile.TemporaryFile(mode="w+") as circuitFile:
  circuitFile.write(charge_amp_test_cap)
  circuitFile.flush()
  circuitFile.seek(0)
  sa = SpiceAnalyzer(circuitFile)
  sa.analyzeManyTrans("Charge_Amplifier_Test_Cap_singletrans.png",1,0,"3",
        [
            "PULSE(0,1,100u,0,0,10u,1)",
            "PULSE(0,1,100u,0,0,100u,1)",
            "PULSE(0,1,100u,0,0,1m,1)",
        ],
        "1u",0,"2m",current=False,debug=False)
