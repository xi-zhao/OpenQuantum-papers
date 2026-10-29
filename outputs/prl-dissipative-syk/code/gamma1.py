import sys, json, numpy as np
sys.path.insert(0, '.')
from pole import pole
from scipy.optimize import brentq
q = 4
mus = np.round(np.concatenate([np.arange(0.0, 0.5, 0.01), np.arange(0.5, 1.21, 0.05)]), 3)
res = []
z_prev = None
for mu in mus:
    z, t, C = pole(q, mu, T=60, dt=0.004, z0=z_prev)
    # also check the slower real branch explicitly (largest Re among roots)
    res.append([float(mu), float(z.real), float(z.imag)])
    z_prev = z
    print(mu, z, flush=True)
json.dump(res, open('../data/gamma1_q4.json', 'w'))
