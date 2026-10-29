import sys, numpy as np
sys.path.insert(0, '.')
import ladder_phys as lp
from ladder2 import Ladder
from pole import pole
from scipy.optimize import brentq
q = int(sys.argv[1])
for mu in [float(x) for x in sys.argv[2].split(',')]:
    z, _, _ = pole(q, mu); thr = -2 * z.real / mu
    ph = lp.Phys.__new__(lp.Phys); ph.Ld = Ladder(q, mu, float(sys.argv[3]), 0.05); ph.n = ph.Ld.n
    es = np.linspace(0.3, 0.999 * thr, 40)
    L = np.array([ph.evals(e * mu, k=6)[:3] for e in es])
    roots = []
    for j in range(3):
        for i in range(len(es) - 1):
            if (L[i, j] - 1) * (L[i + 1, j] - 1) < 0:
                roots.append((j, round(brentq(lambda e: ph.evals(e * mu, k=6)[j] - 1, es[i], es[i + 1], xtol=1e-6), 4)))
    print(f'q={q} mu={mu} 2G1/mu={thr:.3f} crossings {roots}', flush=True)
