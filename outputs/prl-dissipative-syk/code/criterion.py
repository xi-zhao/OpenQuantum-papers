r"""mu->0 rates (units of mu) at finite N: slowest slow (secular, energy-diagonal) and slowest fast
(coherence) modes, from first-order degenerate perturbation theory = eigenvalues of the size
superoperator n restricted to each Bohr-frequency block.  Exact energy mode q, PH mode N-q."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian
q = 4; out = {}
for N in [int(x) for x in sys.argv[1].split(',')]:
    recs = []
    for s in range(int(sys.argv[2])):
        L = EvenLindbladian(N, seed=40 + s); D = L.D
        Ep, Vp = np.linalg.eigh(L.Hp); Em, Vm = np.linalg.eigh(L.Hm)
        C = [Vp.conj().T @ c.toarray() @ Vm for c in L.c]          # <E+|c|E->
        # basis: |a+><b+| (index a,b) and |a-><b-|; n(A)=N/2 A - sum c B c^dag ; n(B)=N/2 B - sum c^dag A c
        # Bohr frequency of |a+><b+| = Ep[a]-Ep[b]; group by frequency (tolerance)
        wP = (Ep[:, None] - Ep[None, :]).ravel(); wM = (Em[:, None] - Em[None, :]).ravel()
        w = np.concatenate([wP, wM]); order = np.argsort(w); ws = w[order]
        cuts = np.where(np.diff(ws) > 1e-7)[0] + 1; groups = np.split(order, cuts)
        D2 = D * D
        def col(idx):
            # n applied to basis element idx, returned as dict over coupled indices restricted later
            v = np.zeros(2 * D2, complex)
            if idx < D2:
                a, b = divmod(idx, D); A = np.zeros((D, D), complex); A[a, b] = 1
                B = -sum(c.conj().T @ A @ c for c in C); v[:D2] = 0.5 * N * A.ravel(); v[D2:] = B.ravel()
            else:
                a, b = divmod(idx - D2, D); B = np.zeros((D, D), complex); B[a, b] = 1
                A = -sum(c @ B @ c.conj().T for c in C); v[D2:] = 0.5 * N * B.ravel(); v[:D2] = A.ravel()
            return v
        slow_min, fast_min = np.inf, np.inf
        for g in groups:
            M = np.array([col(i)[g] for i in g]).T
            ev = np.sort(np.linalg.eigvals(M).real)
            if abs(w[g[0]]) < 1e-7:
                ev = ev[ev > 1e-8]
                # exclude the exact energy mode (rate q) once
                k = np.argmin(np.abs(ev - q)); ev = np.delete(ev, k)
                slow_min = min(slow_min, ev.min())
            else:
                fast_min = min(fast_min, ev.min())
        recs.append([slow_min, fast_min])
    r = np.array(recs)
    out[N] = dict(slow=r[:, 0].tolist(), fast=r[:, 1].tolist())
    print(N, 'second slow (excl. H)/mu:', np.round(r[:, 0], 3), ' slowest fast/mu:', np.round(r[:, 1], 3), ' N/2=', N / 2, flush=True)
json.dump(out, open('../data/criterion.json', 'w'))
