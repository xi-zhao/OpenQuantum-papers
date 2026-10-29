import sys, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian
for N in [8, 10]:
  for s in [40, 41]:
    L = EvenLindbladian(N, seed=s); n = 2 * L.D ** 2
    A = np.empty((n, n), complex); e = np.zeros(n, complex)
    for mu in [0.002, 0.02]:
        for j in range(n):
            e[j] = 1; A[:, j] = L.apply(e, mu); e[j] = 0
        w, V = np.linalg.eig(A); o = np.argsort(-w.real); w = w[o]; V = V[:, o]
        Ep, Vp = np.linalg.eigh(L.Hp); Em, Vm = np.linalg.eigh(L.Hm)
        for k in range(1, 4):
            a, b = L.split(V[:, k]); da = np.diag(Vp.conj().T @ a @ Vp); db = np.diag(Vm.conj().T @ b @ Vm)
            wd = (np.sum(abs(da) ** 2) + np.sum(abs(db) ** 2)) / np.vdot(V[:, k], V[:, k]).real
            print(N, s, mu, k, 'rate/mu=%.4f im=%.2e diagw=%.3f' % (-w[k].real / mu, w[k].imag, wd))
    print(N, s, 'min level spacing +:', np.min(np.diff(Ep)), ' -:', np.min(np.diff(Em)), ' Ep~Em overlap:', np.min(np.abs(Ep[:, None] - Em[None, :])))
