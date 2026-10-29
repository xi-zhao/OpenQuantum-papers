import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from scipy.optimize import brentq
for mu, a, b in [(0.01, 18, 25), (0.1, 6.5, 8.5)]:
    for dt in [0.05, 0.025]:
        ph = Phys(mu, 60, dt)
        f = lambda e: (lambda r: r[np.argmin(np.abs(r - 1))] - 1)(ph.evals(e * mu))
        print(mu, dt, brentq(f, a, b, xtol=1e-6), flush=True)
