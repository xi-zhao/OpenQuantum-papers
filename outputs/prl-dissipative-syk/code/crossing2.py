r"""Finite-N mu_c with mode identification: crossing of q*mu with the slowest FAST mode
(energy-diagonal weight w < 1/2) on span{1,H}^perp."""
import sys, json, time, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian, leading
N = int(sys.argv[1]); seeds = [int(s) for s in sys.argv[2].split(',')]; q = 4
lo0, hi0 = [float(x) for x in sys.argv[3].split(",")]
res = []
for s in seeds:
    L = EvenLindbladian(N, seed=1000 + s)
    Ep, Vp = np.linalg.eigh(L.Hp); Em, Vm = np.linalg.eigh(L.Hm)
    def wdiag(v):
        A, B = L.split(v)
        a = np.diag(Vp.conj().T @ A @ Vp); b = np.diag(Vm.conj().T @ B @ Vm)
        return float((np.sum(abs(a) ** 2) + np.sum(abs(b) ** 2)) / np.vdot(v, v).real)
    log = {}
    def f(mu):
        v, V = leading(L, mu, k=8, tol=1e-7)
        fast = [(-x.real, wdiag(V[:, j])) for j, x in enumerate(v)]
        slowest_any = fast[0]
        fr = [r for r, w in fast if w < 0.5]
        g = min(fr) if fr else np.inf
        log[mu] = dict(fast=g, slowest=slowest_any)
        return g / mu - q
    a, b = lo0, hi0; fa, fb = f(a), f(b); t0 = time.time()
    for _ in range(8):
        m = 0.5 * (a + b); fm = f(m)
        if fm > 0: a, fa = m, fm
        else: b, fb = m, fm
    muc = a - fa * (b - a) / (fb - fa)
    # identity of slowest V_perp mode at mu_c
    sl = log[min(log, key=lambda k: abs(k - muc))]['slowest']
    res.append(dict(N=N, seed=s, mu_c=muc, slowest_at_muc=sl))
    print(N, s, f'mu_c(fast)={muc:.4f}', 'slowest V_perp mode near mu_c: rate/mu=%.3f w=%.2f' % (sl[0] / muc, sl[1]), f'{time.time()-t0:.0f}s', flush=True)
    json.dump(res, open(f'../data/crossing2_N{N}.json', 'w'))
