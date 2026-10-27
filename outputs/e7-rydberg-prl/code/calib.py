"""Light-particle gaps at the working point and at V2 shifted by the reported calibration uncertainty."""
import json
import sys
from pathlib import Path

import numpy as np

from masses import levels
from ryd import PeriodicChain

DC, VC = -1.5011, -1.6260
TH, T = 110.0, 0.2
D0 = DC + T * np.cos(np.deg2rad(TH))
V0 = VC + T * np.sin(np.deg2rad(TH))
SCALES = (-0.017, -0.0085, 0.0, 0.0085, 0.017)

out = []
for L in (int(x) for x in sys.argv[1].split(',')):
  pc = PeriodicChain(L)
  for s in SCALES:
    V = V0 * (1 + s)
    e0, w0, ep, wp = levels(pc, D0, V, n0=6, npi=4)
    out.append(dict(L=L, rel_dV2=s, D=D0, V=V, e0=e0.tolist(), w0=w0.tolist(),
                    ep=ep.tolist(), wp=wp.tolist()))
    print(L, s, ep[0], flush=True)
Path('../data/calib_th110_t0.2.json').write_text(json.dumps(out))
