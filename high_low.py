#!/usr/bin/env python3

import tempfile
from spiceAnalysis import SpiceAnalyzer

circuit = """high_low
r1 1 2 1K
c1 2 0 3n
c2 2 3 3n
r2 3 0 1K
.end
"""

with tempfile.TemporaryFile(mode="w+") as circuitFile:
  circuitFile.write(circuit)
  circuitFile.flush()
  circuitFile.seek(0)
  sa = SpiceAnalyzer(circuitFile)
  sa.analyzeAC("HighLow.png",1,0,["vm(3)"],1,"1k","1000k",debug=False)
  sa.analyzeManyTrans("HighLow2.png",1,0,"vm(3)",
        [
            "PULSE(0,1,10u,0,0,1u,100u)",
            "PULSE(0,1,10u,0,0,10u,100u)",
            "PULSE(0,1,10u,0,0,30u,100u)",
            "PULSE(0,1,10u,0,0,200u,100u)",
        ],
        "0.25u",0,"100u",debug=False)
