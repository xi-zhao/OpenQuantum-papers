import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
q = 4
for mu, es in [(0.1, [14.0, 14.3, 14.6]), (0.15, [9.4, 9.6, 9.9, 10.3]), (0.02, [30, 50, 64])]:
    for S in [20, 40, 60]:
        ph = Phys(mu, S, 0.05)
        print(mu, S, [(e, round(ph.evals(e * mu)[0], 4), round(ph.evals(e * mu)[1], 4)) for e in es], flush=True)
