"""Leading resonance of the large-N fermion two-point function C(t) via the
analytically continued Laplace-domain Dyson equation z + mu + K(z) = 0."""
import sys, numpy as np
sys.path.insert(0, '.')
from sd import solve
from scipy.optimize import fsolve

def pole(q, mu, J=1.0, T=60.0, dt=0.004, z0=None):
    t, C = solve(q=q, J=J, mu=mu, T=T, dt=dt)
    K = 2.0 ** (2 - q) * J**2 * (2 * C) ** (q - 1)
    w = np.ones_like(t) * dt; w[0] = w[-1] = dt / 2
    def F(x):
        z = x[0] + 1j * x[1]
        r = z + mu + np.sum(w * K * np.exp(-z * t))
        return [r.real, r.imag]
    if z0 is None:
        z0 = -0.6 - 0.6 * 1j
    best = None
    for guess in [z0, complex(z0.real, -z0.imag), -mu - 0.5, -(mu + 0.3)]:
        x, info, ier, msg = fsolve(F, [guess.real, guess.imag], full_output=True)
        if ier == 1:
            z = x[0] + 1j * abs(x[1])
            if best is None or z.real > best.real:
                best = z
    return best, t, C

if __name__ == '__main__':
    q = int(sys.argv[1])
    for mu in [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 0.8, 1.2]:
        z, t, C = pole(q, mu)
        print(q, mu, z, 'Gamma1=', -z.real, '2G1/mu=', (-2 * z.real / mu if mu else None), flush=True)
