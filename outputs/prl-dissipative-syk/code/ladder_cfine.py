import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from ladder_complex import spec
from scipy.optimize import minimize
q = 4
mu = float(sys.argv[1]); x0, y0 = float(sys.argv[2]), float(sys.argv[3])
for S in [float(x) for x in sys.argv[4].split(",")]:
    ph = Phys(mu, S, 0.05)
    f = lambda p: np.min(np.abs(spec(ph, p[0] + 1j * p[1]) - 1))
    r = minimize(f, [x0, y0], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-6, maxiter=200))
    print(S, r.x, r.fun, flush=True)
