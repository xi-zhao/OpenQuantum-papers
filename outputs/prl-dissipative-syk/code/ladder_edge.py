import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from pole import pole
from scipy.optimize import brentq
q = 4
for mu in [float(x) for x in sys.argv[1].split(',')]:
    z, _, _ = pole(q, mu); thr = -2 * z.real / mu
    for S in [100, 200, 300]:
        ph = Phys(mu, S, 0.05)
        f = lambda e: (lambda r: r[np.argmin(np.abs(r - 1))] - 1)(ph.evals(e * mu))
        es = np.linspace(0.85 * thr, 0.9999 * thr, 16); fs = [f(e) for e in es]
        rts = [brentq(f, es[i], es[i + 1], xtol=1e-7) for i in range(15) if fs[i] * fs[i + 1] < 0]
        print(f'mu={mu} thr={thr:.5f} S={S} lowest root={min(rts) if rts else None}', flush=True)
