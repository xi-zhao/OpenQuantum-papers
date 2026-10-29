r"""Same as ladder_cbelow but excluding the exact energy mode E=q mu: search |E - q mu| > 0.15 mu,
and also Im E > 0.05 region separately."""
import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from ladder_complex import spec
from pole import pole
from scipy.optimize import minimize
q = 4
for mu in [float(x) for x in sys.argv[1].split(',')]:
    z, _, _ = pole(q, mu); G1, W = -z.real, abs(z.imag)
    for S in [40, 70]:
        ph = Phys(mu, S, 0.05)
        cands = []
        for x in np.linspace(0.05, 0.97, 10) * 2 * G1:
            for y in np.linspace(0.08, 2 * W + 0.4, 8):
                cands.append((np.min(np.abs(spec(ph, x + 1j * y) - 1)), x, y))
        cands.sort()
        best = None
        for d, x, y in cands[:3]:
            f = lambda p: np.min(np.abs(spec(ph, p[0] + 1j * abs(p[1])) - 1))
            r = minimize(f, [x, y], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-8, maxiter=150))
            if best is None or r.fun < best.fun: best = r
        print(f'mu={mu} 2G1={2*G1:.4f} S={S} complex-E best min={best.fun:.2e} at E={best.x[0]:.4f}+{abs(best.x[1]):.4f}i', flush=True)
