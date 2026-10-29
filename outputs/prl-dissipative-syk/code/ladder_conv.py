import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
mu = float(sys.argv[1])
for S, dt in [(20, 0.05), (30, 0.05), (20, 0.025), (40, 0.05)]:
    ph = Phys(mu, S, dt)
    print(S, dt, [np.round(ph.evals(e * mu)[:3], 5).tolist() for e in [4.0, 14.0, 14.3, 14.6]], flush=True)
