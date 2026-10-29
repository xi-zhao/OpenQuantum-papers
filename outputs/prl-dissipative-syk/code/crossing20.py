import sys, json, time, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian, leading
N = 20; s = int(sys.argv[1]); q = 4
L = EvenLindbladian(N, q=q, seed=3000 + s)
pts = {}
def f(mu):
    v, V = leading(L, mu, k=2, tol=1e-6)
    pts[mu] = -v[0].real
    print(s, mu, -v[0].real, -v[0].real / mu, flush=True)
    return -v[0].real / mu - q
a, b = 0.18, 0.26
fa, fb = f(a), f(b)
for _ in range(3):
    m = a - fa * (b - a) / (fb - fa); fm = f(m)
    if fm > 0: a, fa = m, fm
    else: b, fb = m, fm
muc = a - fa * (b - a) / (fb - fa)
print('mu_c', muc)
json.dump(dict(N=N, seed=s, mu_c=muc, pts=pts), open(f'../data/crossing20_s{s}.json', 'w'))
