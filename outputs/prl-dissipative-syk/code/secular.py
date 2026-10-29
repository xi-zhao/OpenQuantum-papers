"""Secular (mu->0) slow generator: Pauli master equation on energy eigenstates (even sector),
W_{EE'} = sum_i |<E|psi_i|E'>|^2 - (N/2) delta_{EE'}  (in units of mu).  N = 6 mod 8 or 0 mod 8
to avoid intra-sector degeneracies.  Eigenvalues/q -> k (OU tower) as N -> infinity."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian
out = {}
for N in [int(x) for x in sys.argv[1].split(',')]:
    lams = []
    for s in range(3):
        L = EvenLindbladian(N, seed=9000 + s)
        Wtot = []
        # rho even: blocks (+,+) and (-,-); psi maps + <-> -, so W couples E_+ with E_-
        Ep, Vp = np.linalg.eigh(L.Hp); Em, Vm = np.linalg.eigh(L.Hm)
        T = np.zeros((L.D, L.D))
        for c in L.c:
            M = Vp.conj().T @ c.toarray() @ Vm
            T += np.abs(M) ** 2
        W = np.block([[np.zeros((L.D, L.D)), T], [T.T, np.zeros((L.D, L.D))]]) - 0.5 * N * np.eye(2 * L.D)
        w = np.sort(np.linalg.eigvalsh(W))[::-1]
        lams.append((-w[:7] / 4).tolist())
    lam = np.mean(lams, 0)
    out[N] = lam.tolist()
    print(N, np.round(lam, 3), flush=True)
json.dump(out, open('../data/secular.json', 'w'))
