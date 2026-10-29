r"""Complex bound states below the EP: minimize min_k|lambda_k(E)-1| over complex E with
Re E < 2 Gamma1, |Im E| up to 2|Im z1|+0.3; check stability in S."""
import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from ladder_complex import spec
from pole import pole
from scipy.optimize import minimize
q = 4
for mu in [float(x) for x in sys.argv[1].split(',')]:
    z, _, _ = pole(q, mu); G1, W = -z.real, abs(z.imag)
    res = {}
    for S in [40, 70]:
        ph = Phys(mu, S, 0.05)
        best = (9, 0, 0)
        for x in np.linspace(0.1, 0.97, 9) * 2 * G1:
            for y in np.linspace(0.05, 2 * W + 0.3, 8):
                d = np.min(np.abs(spec(ph, x + 1j * y) - 1))
                if d < best[0]: best = (d, x, y)
        f = lambda p: np.min(np.abs(spec(ph, p[0] + 1j * p[1]) - 1))
        r = minimize(f, [best[1], best[2]], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-7, maxiter=150))
        res[S] = (r.fun, r.x)
        print(f'mu={mu} 2G1={2*G1:.4f} 2W={2*W:.4f} S={S} min={r.fun:.2e} at E={r.x[0]:.4f}+{r.x[1]:.4f}i', flush=True)
