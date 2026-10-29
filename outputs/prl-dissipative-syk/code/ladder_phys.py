r"""Physical (exchange-antisymmetric, X=+1) bound states of the ladder below the 2*Gamma1 threshold.
Bound state of decay rate E <=> K(E) has eigenvalue 1 on the X=+1 subspace."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder2 import Ladder, PAIRS
from pole import pole
from scipy.sparse.linalg import LinearOperator, eigs
from scipy.optimize import brentq
q = 4

class Phys:
    def __init__(self, mu, S=20, dt=0.05):
        self.Ld = Ladder(q, mu, S, dt); self.n = self.Ld.n
    def X(self, v):
        n = self.n; out = np.zeros_like(v)
        for i, (a, b) in enumerate(PAIRS):
            j = PAIRS.index((b, a)); out[i*n:(i+1)*n] = -v[j*n:(j+1)*n][::-1]
        return out
    def evals(self, E, k=6):
        self.Ld.setE(E)
        P = lambda v: 0.5 * (v + self.X(v))
        N4 = 4 * self.n
        op = LinearOperator((N4, N4), matvec=lambda v: P(self.Ld.apply(P(v))), dtype=complex)
        v = eigs(op, k=k, which='LR', tol=1e-9, ncv=40, return_eigenvectors=False)
        return np.sort(v.real[np.abs(v.imag) < 1e-6])[::-1]

def scan(mu, nE=36, S=40, dt=0.05):
    z, _, _ = pole(q, mu); G1 = -z.real
    ph = Phys(mu, S, dt)
    Es = np.linspace(0.02 * mu, 0.995 * 2 * G1, nE)
    tops = [ph.evals(E) for E in Es]
    cnt = [int(np.sum(t > 1)) for t in tops]
    def f(E):
        r = ph.evals(E); return r[np.argmin(np.abs(r - 1))] - 1
    roots = []
    for i in range(nE - 1):
        if cnt[i] != cnt[i + 1]:
            try: roots.append(brentq(f, Es[i], Es[i + 1], xtol=1e-6))
            except ValueError: roots.append(float('nan'))
    return dict(mu=mu, G1=G1, roots=roots, n0=cnt[0], Es=Es.tolist(), tops=[t[:4].tolist() for t in tops])

if __name__ == '__main__':
    out = []
    for mu in [float(x) for x in sys.argv[1].split(',')]:
        r = scan(mu); out.append(r)
        print(f"mu={mu:.3f} 2G1/mu={2*r['G1']/mu:7.3f} n(lam>1,E->0)={r['n0']} bound E/mu =",
              [round(x / mu, 4) for x in r['roots']], flush=True)
        json.dump(out, open(f'../data/ladder_phys_{sys.argv[2]}.json', 'w'))
