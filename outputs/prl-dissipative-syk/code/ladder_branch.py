r"""Correct bound-state search: track sorted physical-sector eigenvalues lambda_j(E), j=1..3,
and find E where each branch crosses 1 (continuous in E).  Bound state <=> lambda_j(E)=1."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from pole import pole
from scipy.optimize import brentq
q = 4
out = []
for mu in [float(x) for x in sys.argv[1].split(',')]:
    z, _, _ = pole(q, mu); thr = -2 * z.real / mu
    S = float(sys.argv[2]) if len(sys.argv) > 2 else 60
    ph = Phys(mu, S, 0.05)
    lam = lambda e, j: ph.evals(e * mu, k=6)[j]
    es = np.linspace(0.3, 0.999 * thr, 50)
    L = np.array([ph.evals(e * mu, k=6)[:3] for e in es])
    roots = []
    for j in range(3):
        for i in range(len(es) - 1):
            if (L[i, j] - 1) * (L[i + 1, j] - 1) < 0:
                roots.append((j, brentq(lambda e: lam(e, j) - 1, es[i], es[i + 1], xtol=1e-6)))
    out.append(dict(mu=mu, thr=thr, S=S, roots=roots, grid=es.tolist(), lam=L.tolist()))
    print(f'mu={mu:.3f} 2G1/mu={thr:.3f} S={S} branch crossings (j, E/mu): {[(j, round(r,4)) for j,r in roots]}', flush=True)
    json.dump(out, open(f'../data/ladder_branch_{sys.argv[3] if len(sys.argv)>3 else "x"}.json', 'w'))
