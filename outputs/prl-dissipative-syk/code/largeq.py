import sys, json, numpy as np
sys.path.insert(0, '.')
from sd import solve
from scipy.optimize import fsolve, brentq

def gamma1(q, mu, T, dt):
    t, C = solve(q=q, mu=mu, T=T, dt=dt)
    K = 2.0 ** (2 - q) * (2 * C) ** (q - 1)
    w = np.full_like(t, dt); w[0] = w[-1] = dt / 2
    def F(x):
        z = x[0] + 1j * x[1]; r = z + mu + np.sum(w * K * np.exp(-z * t)); return [r.real, r.imag]
    Jc = np.sqrt(q / 2 ** (q - 1)); a = np.sqrt(Jc**2 + (q * mu / 2) ** 2)
    best = None
    for g in [-2 * a / q, -2 * a / q - 0.3j * a / q, -mu]:
        x, info, ier, msg = fsolve(F, [g.real, g.imag], full_output=True)
        if ier == 1 and x[0] < 0:
            if best is None or x[0] > best: best = x[0]
    return -best

out = {}
for q in [4, 6, 8, 12]:
    Jc = np.sqrt(q / 2 ** (q - 1))
    pred = 4 * Jc / (q * np.sqrt(q * q - 4))
    T = max(60, 40 * q / Jc); dt = min(0.004 * q / 4 / Jc * 0.7, 0.02)
    f = lambda m: q * m - 2 * gamma1(q, m, T, dt)
    muc = brentq(f, 0.3 * pred, 3 * pred, xtol=1e-5)
    out[q] = [muc, pred]
    print(q, 'mu_c SD', round(muc, 5), 'large-q formula', round(pred, 5), 'ratio', round(muc / pred, 3), flush=True)
json.dump(out, open('../data/largeq_muc.json', 'w'))
