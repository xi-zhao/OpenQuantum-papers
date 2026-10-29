r"""Large-N real-time Schwinger-Dyson (Kadanoff-Baym) equation for the
infinite-temperature steady state of dissipative Majorana SYK.
C(t) = <psi(t) psi(0)>_{beta=0},  C(0)=1/2:
  dC/dt = -mu C - 2^{2-q} J^2 \int_0^t ds (2C(t-s))^{q-1} C(s)
(normalization fixed so that the large-q limit reproduces g'' = -2 Jc^2 e^g,
 g'(0) = -q mu, Jc^2 = q J^2 / 2^{q-1}).  Solved with a 2nd-order predictor-corrector."""
import numpy as np

def solve(q=4, J=1.0, mu=0.0, T=40.0, dt=0.005):
    n = int(T / dt) + 1
    C = np.zeros(n); C[0] = 0.5
    K = np.zeros(n)  # kernel k(t) = 2^{2-q} J^2 (2C)^{q-1}
    pref = 2.0 ** (2 - q) * J ** 2
    K[0] = pref * (2 * C[0]) ** (q - 1)
    def rhs(i):
        # trapezoid of K(t_i - s) C(s) on [0, t_i]
        if i == 0:
            return -mu * C[0]
        w = K[i::-1] * C[: i + 1]
        integ = dt * (w.sum() - 0.5 * (w[0] + w[-1]))
        return -mu * C[i] - integ
    f = np.zeros(n); f[0] = rhs(0)
    for i in range(1, n):
        C[i] = C[i - 1] + dt * f[i - 1]  # predictor
        for _ in range(2):
            K[i] = pref * (2 * C[i]) ** (q - 1)
            f[i] = rhs(i)
            C[i] = C[i - 1] + 0.5 * dt * (f[i - 1] + f[i])
        K[i] = pref * (2 * C[i]) ** (q - 1); f[i] = rhs(i)
    t = np.arange(n) * dt
    return t, C

def rate(t, C, window=(0.5, 0.9)):
    """Asymptotic exponential decay rate of |C| from log-envelope fit."""
    a, b = int(window[0] * len(t)), int(window[1] * len(t))
    y = np.log(np.abs(C[a:b]) + 1e-300)
    # use local maxima of |C| to be robust to oscillations
    idx = [i for i in range(1, len(y) - 1) if y[i] >= y[i - 1] and y[i] >= y[i + 1]]
    if len(idx) < 5:
        idx = list(range(len(y)))
    p = np.polyfit(t[a:b][idx], y[idx], 1)
    return -p[0]
