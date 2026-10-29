import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from pole import pole
from scipy.optimize import brentq
q = 4
mu = float(sys.argv[1]); j = int(sys.argv[2])
z, _, _ = pole(q, mu); thr = -2 * z.real / mu
for S in [float(x) for x in sys.argv[3].split(',')]:
    ph = Phys(mu, S, 0.05)
    es = np.linspace(0.9 * thr if j == 1 else 0.8 * min(4, thr), 0.9999 * thr, 30)
    L = np.array([ph.evals(e * mu, k=6)[j] for e in es])
    rts = [brentq(lambda e: ph.evals(e * mu, k=6)[j] - 1, es[i], es[i + 1], xtol=1e-7)
           for i in range(len(es) - 1) if (L[i] - 1) * (L[i + 1] - 1) < 0]
    print(f'mu={mu} branch {j} thr={thr:.5f} S={S:.0f} roots={np.round(rts, 5).tolist()}  gap-to-thr={[round(thr - r, 5) for r in rts]}', flush=True)
