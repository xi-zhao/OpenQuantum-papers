r"""Collective bound state E_b(mu) of the physical singlet ladder below 2 Gamma1,
excluding the exact energy mode at q mu.  Roots of lambda(E)=1 in (lo, 0.999*2G1) for S = 40,60,80,
then S->inf by fitting E(S) = E_inf + a e^{-S/xi} or E_inf - c/S (reported both)."""
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
    for S in [40, 60, 80, 100]:
        ph = Phys(mu, S, 0.05)
        f = lambda e: (lambda r: r[np.argmin(np.abs(r - 1))] - 1)(ph.evals(e * mu))
        lo = max(4.3, 0.6 * thr) if thr > 4.6 else 0.6 * thr
        es = np.linspace(lo, 0.9995 * thr, 18); fs = [f(e) for e in es]
        rts = [brentq(f, es[i], es[i + 1], xtol=1e-6) for i in range(17) if fs[i] * fs[i + 1] < 0]
        rts = [r for r in rts if abs(r - 4) > 0.02]
        rows[S] = rts
        print(f'mu={mu:.3f} thr={thr:.4f} S={S} roots={np.round(rts, 4).tolist()}', flush=True)
    out.append(dict(mu=mu, thr=thr, roots={str(k): v for k, v in rows.items()}))
    json.dump(out, open(f'../data/ladder_Eb_{sys.argv[2]}.json', 'w'))
