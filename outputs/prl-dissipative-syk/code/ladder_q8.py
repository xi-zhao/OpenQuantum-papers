import sys, numpy as np
sys.path.insert(0, '.')
import ladder_phys as lp
from ladder2 import Ladder
from pole import pole
from scipy.optimize import brentq
q = int(sys.argv[1]); mu = float(sys.argv[2])
z, _, _ = pole(q, mu); thr = -2 * z.real / mu
for S in [float(x) for x in sys.argv[3].split(',')]:
    ph = lp.Phys.__new__(lp.Phys); ph.Ld = Ladder(q, mu, S, 0.05); ph.n = ph.Ld.n
    es = np.linspace(0.5 * min(q, thr), 0.9995 * thr, 30)
    L = np.array([ph.evals(e * mu, k=6)[0] for e in es])
    rts = [brentq(lambda e: ph.evals(e * mu, k=6)[0] - 1, es[i], es[i + 1], xtol=1e-6) for i in range(29) if (L[i] - 1) * (L[i + 1] - 1) < 0]
    print(f'q={q} mu={mu} thr={thr:.4f} S={S:.0f} roots={np.round(rts,4).tolist()}', flush=True)
