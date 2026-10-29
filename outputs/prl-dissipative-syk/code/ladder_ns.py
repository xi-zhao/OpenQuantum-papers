r"""Non-singlet channel: K_ns = -K/(q-1); bound state <=> K has eigenvalue -(q-1) on X=+1 subspace."""
import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from pole import pole
from scipy.sparse.linalg import LinearOperator, eigs
q = 4
for mu in [float(x) for x in sys.argv[1].split(',')]:
    z, _, _ = pole(q, mu); G1 = -z.real
    ph = Phys(mu)
    mins = []
    for E in np.linspace(0.02 * mu, 0.99 * 2 * G1, 12):
        ph.Ld.setE(E)
        P = lambda v: 0.5 * (v + ph.X(v)); N4 = 4 * ph.n
        op = LinearOperator((N4, N4), matvec=lambda v: P(ph.Ld.apply(P(v))), dtype=complex)
        v = eigs(op, k=4, which='SR', tol=1e-8, ncv=40, return_eigenvectors=False)
        mins.append(v.real.min())
    print(f'mu={mu} min Re eig(K) over E in (0,2G1): {min(mins):.4f}  (bound state needs {-(q-1)})', flush=True)
