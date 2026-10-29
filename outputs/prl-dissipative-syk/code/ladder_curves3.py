import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from pole import pole
q = 4; out = {}
for mu in [0.1, 0.25, 0.4]:
    z, _, _ = pole(q, mu); thr = -2 * z.real
    ph = Phys(mu, 100, 0.05)
    Es = np.linspace(0.02, 0.995, 34) * thr
    lam = [ph.evals(E, k=6)[:3].tolist() for E in Es]
    out[mu] = dict(thr=thr, E=Es.tolist(), lam=lam)
    print(mu, 'done', flush=True)
json.dump(out, open('../data/ladder_curves3.json', 'w'))
