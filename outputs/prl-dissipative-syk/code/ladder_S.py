import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from scipy.optimize import brentq
q = 4
for mu in [float(x) for x in sys.argv[1].split(',')]:
    for S in [20, 40, 60, 80]:
        ph = Phys(mu, S, 0.05)
        f = lambda E: (lambda r: r[np.argmin(np.abs(r - 1))] - 1)(ph.evals(E))
        try:
            r = brentq(f, 3.0 * mu, 4.3 * mu, xtol=1e-7)
        except ValueError as e:
            r = float('nan')
        print(f'mu={mu} S={S} E_b/mu={r/mu:.4f}', flush=True)
