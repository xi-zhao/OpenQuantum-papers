"""mu_c(N): where the slowest energy-orthogonal even mode crosses the exact energy mode q*mu."""
import sys, json, time, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian, leading
N = int(sys.argv[1]); seeds = [int(s) for s in sys.argv[2].split(',')]
q = 4
res = []
for s in seeds:
    L = EvenLindbladian(N, q=q, seed=1000 + s)
    cache = {}
    def f(mu):
        v, V = leading(L, mu, k=3, tol=1e-7)
        cache[mu] = (v, L.mean_size(V[:, 0]))
        return -v[0].real / mu - q   # >0 : energy mode is the gap
    a, b = 0.06, 0.6
    fa, fb = f(a), f(b)
    t0 = time.time()
    for _ in range(9):
        m = 0.5 * (a + b); fm = f(m)
        if fm > 0: a, fa = m, fm
        else: b, fb = m, fm
    muc = a - fa * (b - a) / (fb - fa)
    rec = dict(N=N, seed=s, mu_c=muc, scan={str(k): [-v[0].real, v[0].imag, sz] for k, (v, sz) in cache.items()})
    res.append(rec)
    print(N, s, 'mu_c=%.4f' % muc, '%.0fs' % (time.time() - t0), flush=True)
    json.dump(res, open(f'../data/crossing_N{N}.json', 'w'))
