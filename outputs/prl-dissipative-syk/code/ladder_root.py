import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder import greens, kernel
from scipy.linalg import eigvals
from scipy.optimize import brentq
q = 4
def top(mu, E, Gs, m=4):
    s, K = kernel(q, mu, E, G=Gs)
    v = eigvals(K)
    v = v[np.abs(v.imag) < 1e-6].real
    return np.sort(v)[::-1][:m]
mu = float(sys.argv[1])
for S, dt in [(20, 0.1), (25, 0.07), (30, 0.05)]:
    Gs = greens(q, mu, S, dt)
    f = lambda E: top(mu, E, Gs, 1)[0] - 1
    Es = np.linspace(0.5 * q * mu, 1.4 * q * mu, 10)
    vals = [f(E) for E in Es]
    roots = [brentq(f, Es[i], Es[i + 1], xtol=1e-5) for i in range(len(Es) - 1) if vals[i] * vals[i + 1] < 0]
    print(S, dt, 'roots E/(q mu)=', [round(r / (q * mu), 4) for r in roots], flush=True)
