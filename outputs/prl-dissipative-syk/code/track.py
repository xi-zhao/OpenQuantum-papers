r"""ED with mode identification.  Even sector; for each mu compute k leading eigenpairs of L on
span{1,H}^perp and the energy-diagonal weight w = ||D(X)||^2/||X||^2 (D: dephasing in the H eigenbasis,
per parity block).  Slow modes: w > 1/2; fast modes: w < 1/2.  Reports slowest slow and slowest fast."""
import sys, json, time, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian, leading
N = int(sys.argv[1]); seeds = [int(s) for s in sys.argv[2].split(',')]
mus = [float(x) for x in sys.argv[3].split(',')]
k = int(sys.argv[4]) if len(sys.argv) > 4 else 6
out = []
for s in seeds:
    L = EvenLindbladian(N, seed=5000 + s)
    Ep, Vp = np.linalg.eigh(L.Hp); Em, Vm = np.linalg.eigh(L.Hm)
    def wdiag(v):
        A, B = L.split(v)
        a = np.diag(Vp.conj().T @ A @ Vp); b = np.diag(Vm.conj().T @ B @ Vm)
        return float((np.sum(np.abs(a) ** 2) + np.sum(np.abs(b) ** 2)) / np.vdot(v, v).real)
    for mu in mus:
        t = time.time()
        v, V = leading(L, mu, k=k, tol=1e-7)
        recs = [[-x.real, x.imag, wdiag(V[:, j]), L.mean_size(V[:, j])] for j, x in enumerate(v)]
        slow = [r for r in recs if r[2] > 0.5]; fast = [r for r in recs if r[2] <= 0.5]
        out.append(dict(seed=s, mu=mu, modes=recs))
        print(N, s, mu, 'slow', [round(r[0] / mu, 3) for r in slow][:3], 'fast', [round(r[0] / mu, 3) for r in fast][:2],
              'w', [round(r[2], 2) for r in recs], f'{time.time()-t:.0f}s', flush=True)
        json.dump(out, open(f'../data/track_N{N}.json', 'w'))
