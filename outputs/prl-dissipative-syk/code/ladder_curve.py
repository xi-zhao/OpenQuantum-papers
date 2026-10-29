r"""Top physical-sector (X=+1) ladder eigenvalues lambda(E) on a fine E grid; also the extent
<s^2>^{1/2} of the eigenvector closest to 1 (bound-state size in relative time)."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from scipy.sparse.linalg import LinearOperator, eigs
mu = float(sys.argv[1]); S = float(sys.argv[2]); emax = float(sys.argv[3])
ph = Phys(mu, S, 0.05); n = ph.n; s = ph.Ld.s
out = []
for e in np.linspace(0.5, emax, int(sys.argv[4]) if len(sys.argv) > 4 else 30):
    ph.Ld.setE(e * mu)
    P = lambda v: 0.5 * (v + ph.X(v))
    op = LinearOperator((4 * n,) * 2, matvec=lambda v: P(ph.Ld.apply(P(v))), dtype=complex)
    w, V = eigs(op, k=4, which='LR', tol=1e-9, ncv=40)
    o = np.argsort(-w.real); w = w[o]; V = V[:, o]
    j = np.argmin(np.abs(w - 1)); v = V[:, j]
    dens = sum(np.abs(v[i * n:(i + 1) * n]) ** 2 for i in range(4))
    ext = np.sqrt((dens * s ** 2).sum() / dens.sum())
    out.append([e, *w.real[:3].tolist(), float(ext)])
    print(f'E/mu={e:6.2f}  lam={np.round(w.real[:3],4)}  closest-to-1 extent={ext:.2f}', flush=True)
json.dump(out, open(f'../data/ladder_curve_mu{mu}.json', 'w'))
