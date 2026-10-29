r"""Search for complex-E bound states (oscillatory resonances) of the physical singlet ladder
with Re E < 2 Gamma1: min_k |lambda_k(E) - 1| on a grid; also the non-singlet channel
(lambda = -(q-1))."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from pole import pole
from scipy.sparse.linalg import LinearOperator, eigs
q = 4
def spec(ph, E, k=10):
    ph.Ld.setE(E)
    P = lambda v: 0.5 * (v + ph.X(v)); N4 = 4 * ph.n
    op = LinearOperator((N4, N4), matvec=lambda v: P(ph.Ld.apply(P(v))), dtype=complex)
    return eigs(op, k=k, which='LM', tol=1e-7, ncv=50, return_eigenvectors=False)
if __name__ == "__main__":
    out = []
    for mu in [float(x) for x in sys.argv[1].split(',')]:
        z, _, _ = pole(q, mu); G1 = -z.real
        ph = Phys(mu, 30, 0.05)
        best = (9, None, None)
        for x in np.linspace(0.05, 0.97, 10) * 2 * G1:
            for y in np.linspace(0.1, 1.6, 9):
                v = spec(ph, x + 1j * y)
                d1 = np.min(np.abs(v - 1)); d2 = np.min(np.abs(v + (q - 1)))
                d = min(d1, d2)
                if d < best[0]: best = (d, x, y)
        out.append(dict(mu=mu, mind=best[0], at=[best[1], best[2]]))
        print(f'mu={mu} min|lam-target| over complex grid = {best[0]:.3f} at E=({best[1]:.3f},{best[2]:.3f})', flush=True)
    json.dump(out, open('../data/ladder_complex.json', 'w'))
