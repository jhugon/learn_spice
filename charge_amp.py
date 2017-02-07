#!/usr/bin/env python3

import tempfile
from spiceAnalysis import SpiceAnalyzer

library = """*
* 1 input+, 2 input-, 3 output, 4 ground
.subckt opamp 1 2 3 4
E1 3 4 1 2 1e8
.ends
*
*
*
* 1e-4 wT, 1e-3 wF, 1e5 gain
* 1 input+, 2 input-, 3 output, 4 ground
.subckt charge_amp 2 0 5 0
rf 2 3 100k
cf 2 3 10n
xin 0 2 3 0 opamp
rc 3 4 1Meg
cc 3 4 1n
rt 4 5 10k
ct 4 5 10n
xt 0 4 5 0 opamp
.ends
"""

charge_amp = """charge_amp
*""" + library + """*
*
ri 1 0 1Meg
ci 1 0 100p
* blocking cap or not
cb 1 2 1u
*rb 1 2 0
*
xca 2 0 3 0 charge_amp
.end
"""

charge_amp_test_cap = """charge_amp_test_cap
*""" + library + """*
*
ct 1 2 100n
rt 1 2 1000G
*
xca 2 0 3 0 charge_amp
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
              "PULSE(0,1n,0,0,0,1u,1)",
              "PULSE(0,1n,0,0,0,2u,1)",
              "PULSE(0,2n,0,0,0,2u,1)",
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
  sa.analyzeManyTrans(savename+"Charge_Amplifier_Test_Cap_singletrans.png",1,0,"3",
        [
            "PULSE(0,1,0,0,0,10u,1)",
            "PULSE(0,1,0,0,0,100u,1)",
            "PULSE(0,1,0,0,0,1m,1)",
        ],
        "1u",0,"2m",current=False,debug=False)
