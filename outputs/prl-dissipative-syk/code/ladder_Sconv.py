import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from pole import pole
from scipy.optimize import brentq
q = 4
out = {}
for mu in [float(x) for x in sys.argv[1].split(',')]:
    z, _, _ = pole(q, mu); thr = -2 * z.real / mu
    row = []
    for S in [20, 40, 60, 80]:
        ph = Phys(mu, S, 0.05)
        f = lambda e: (lambda r: r[np.argmin(np.abs(r - 1))] - 1)(ph.evals(e * mu))
        lo, hi = 0.75 * min(4, thr), 0.999 * thr
        es = np.linspace(lo, hi, 14); fs = [f(e) for e in es]
        rts = [brentq(f, es[i], es[i+1], xtol=1e-6) for i in range(13) if fs[i] * fs[i+1] < 0]
        row.append(rts)
        print(f'mu={mu} thr(2G1/mu)={thr:.4f} S={S} roots E/mu={np.round(rts,4).tolist()}', flush=True)
    out[mu] = dict(thr=thr, roots=row)
json.dump(out, open(f'../data/ladder_Sconv_{sys.argv[2]}.json', 'w'))
