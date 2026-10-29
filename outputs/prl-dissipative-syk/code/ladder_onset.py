import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from pole import pole
from scipy.optimize import brentq
q = 4
out = []
for mu in [float(x) for x in sys.argv[1].split(',')]:
    z, _, _ = pole(q, mu); thr = -2 * z.real / mu
    rows = {}
    for S, dt in [(60, 0.05), (100, 0.05)]:
        ph = Phys(mu, S, dt)
        f = lambda e: (lambda r: r[np.argmin(np.abs(r - 1))] - 1)(ph.evals(e * mu))
        es = np.linspace(2.05, 0.9995 * thr, 60); fs = [f(e) for e in es]
        rts = [brentq(f, es[i], es[i + 1], xtol=1e-6) for i in range(21) if fs[i] * fs[i + 1] < 0]
        rows[f'{S}_{dt}'] = rts
        print(f'mu={mu:.3f} thr={thr:.4f} S={S} dt={dt} roots={np.round(rts, 4).tolist()}', flush=True)
    out.append(dict(mu=mu, thr=thr, roots=rows))
    json.dump(out, open(f'../data/ladder_onset_{sys.argv[2]}.json', 'w'))
