r"""mu->0 at finite N: rates of modes outside ker[H,.] should be mu*N/2 (+O(mu^2)) when Bohr frequencies
are non-degenerate.  Dense diagonalization; slow sector = ker[H,.] via projection onto energy-degenerate blocks."""
import sys, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian
for N in [8, 10, 12]:
    L = EvenLindbladian(N, seed=40); n = 2 * L.D ** 2
    Ep, Vp = np.linalg.eigh(L.Hp); Em, Vm = np.linalg.eigh(L.Hm)
    def kerw(v, tol=1e-8):
        A, B = L.split(v)
        a = Vp.conj().T @ A @ Vp; b = Vm.conj().T @ B @ Vm
        ma = np.abs(Ep[:, None] - Ep[None, :]) < tol; mb = np.abs(Em[:, None] - Em[None, :]) < tol
        return (np.sum(np.abs(a[ma]) ** 2) + np.sum(np.abs(b[mb]) ** 2)) / np.vdot(v, v).real
    for mu in [0.001, 0.01]:
        A = np.empty((n, n), complex); e = np.zeros(n, complex)
        for j in range(n):
            e[j] = 1; A[:, j] = L.apply(e, mu); e[j] = 0
        w, V = np.linalg.eig(A)
        ws = np.array([kerw(V[:, j]) for j in range(n)])
        r = -w.real / mu
        fast = r[ws < 0.5]; slow = r[(ws >= 0.5) & (r > 1e-6)]
        print(f'N={N} mu={mu} fast: min={fast.min():.3f} median={np.median(fast):.3f} (N/2={N/2})  slow: min={slow.min():.3f}', flush=True)
