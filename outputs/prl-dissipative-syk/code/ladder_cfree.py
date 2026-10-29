r"""Complex bound-state search off the real axis: minimize min_k|lambda_k(E)-1| with Im E >= ymin."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from ladder_complex import spec
from pole import pole
from scipy.optimize import minimize
q = 4; ymin = 0.05
out = []
for mu in [float(x) for x in sys.argv[1].split(',')]:
    z, _, _ = pole(q, mu); G1, W = -z.real, abs(z.imag)
    for S in [40, 70]:
        ph = Phys(mu, S, 0.05)
        cands = []
        for x in np.linspace(0.05, 0.98, 12) * 2 * G1:
            for y in np.linspace(ymin, 2 * W + 0.6, 9):
                cands.append((np.min(np.abs(spec(ph, x + 1j * y) - 1)), x, y))
        cands.sort(); best = None
        for d, x, y in cands[:4]:
            f = lambda p: np.min(np.abs(spec(ph, p[0] + 1j * max(p[1], ymin)) - 1)) + 10 * max(0, ymin - p[1])
            r = minimize(f, [x, y], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-8, maxiter=120))
            if best is None or r.fun < best.fun: best = r
        out.append(dict(mu=mu, S=S, min=best.fun, E=[best.x[0], best.x[1]], twoG1=2 * G1))
        print(f'mu={mu} 2G1={2*G1:.4f} S={S} min|lam-1| (ImE>={ymin}) = {best.fun:.2e} at E={best.x[0]:.4f}+{best.x[1]:.4f}i', flush=True)
        json.dump(out, open('../data/ladder_cfree.json', 'w'))
