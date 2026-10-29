r"""Matrix-free ladder kernel for the vectorized dissipative SYK (two-site picture).
Psi_ab(s), s = t1 - t2 in [-S, S]; channel factor 1 (singlet) or -1/(q-1) (non-singlet).
(K phi)_ab(s) = sum_cd int dr M_{ab,cd}(s - r) w_cd(r) phi_cd(r),
M_{ab,cd}(x) = int du g_ac(u) g_bd(u - x),  g = G e^{E u/2},
w_cd(r) = -(q-1) J^2 c_c c_d G_cd(r)^{q-2}."""
import sys, numpy as np
sys.path.insert(0, '.')
from sd import solve
from scipy.signal import fftconvolve
from scipy.sparse.linalg import LinearOperator, eigs

FL = ['L', 'R']; PAIRS = [(a, b) for a in FL for b in FL]
cf = {'L': -1j, 'R': 1j}

class Ladder:
    def __init__(self, q, mu, S, dt, J=1.0):
        self.q, self.mu, self.S, self.dt = q, mu, S, dt
        t, C = solve(q=q, J=J, mu=mu, T=S + 2 * dt, dt=min(dt, 0.004))
        n = int(round(S / dt)); self.n = 2 * n + 1
        s = np.arange(-n, n + 1) * dt; self.s = s
        Ca = np.interp(np.abs(s), t, C); sg = np.sign(s)
        self.G = {('L', 'L'): sg * Ca, ('R', 'R'): sg * Ca, ('L', 'R'): 1j * Ca, ('R', 'L'): -1j * Ca}
        # G^{q-2}: same-side sign^{q-2}=1 including s=0
        Gp = {('L', 'L'): Ca ** (q - 2), ('R', 'R'): Ca ** (q - 2),
              ('L', 'R'): (1j * Ca) ** (q - 2), ('R', 'L'): (-1j * Ca) ** (q - 2)}
        self.w = {cd: -(q - 1) * J ** 2 * cf[cd[0]] * cf[cd[1]] * Gp[cd] for cd in PAIRS}

    def setE(self, E, channel=1.0):
        dt, n = self.dt, self.n
        ew = np.exp(E * self.s / 2)
        g = {k: v * ew for k, v in self.G.items()}
        self.M = {}
        for (a, b) in PAIRS:
            for (c, d) in PAIRS:
                self.M[(a, b, c, d)] = fftconvolve(g[(a, c)], g[(b, d)][::-1], mode='full') * dt * channel
        # pad length for linear convolution with phi
        self.E = E

    def apply(self, v):
        n, dt = self.n, self.dt
        phi = {cd: v[i * n:(i + 1) * n] for i, cd in enumerate(PAIRS)}
        out = []
        for (a, b) in PAIRS:
            acc = np.zeros(n, complex)
            for (c, d) in PAIRS:
                f = self.w[(c, d)] * phi[(c, d)]
                conv = fftconvolve(self.M[(a, b, c, d)], f, mode='full') * dt
                # M index k <-> x=(k-(n-1))dt ; result index m <-> s=(m-(n-1)-(n-1))dt + s_min...
                # s_i = (i - (n-1)/2)dt ; r_j likewise; s_i - r_j = (i-j)dt -> M index i-j+n-1
                # conv[m] = sum_j M[m-j] f[j]  -> want m-j = i-j+n-1 -> m = i+n-1
                acc += conv[n - 1:2 * n - 1]
            out.append(acc)
        return np.concatenate(out)

    def eig(self, k=6, which='LM', sigma=None):
        N4 = 4 * self.n
        op = LinearOperator((N4, N4), matvec=self.apply, dtype=complex)
        vals, vecs = eigs(op, k=k, which=which, tol=1e-9, ncv=max(3 * k, 30))
        return vals, vecs

if __name__ == '__main__':
    q = 4; mu = float(sys.argv[1])
    for S, dt in [(20, 0.1), (20, 0.05), (20, 0.025)]:
        Ld = Ladder(q, mu, S, dt); Ld.setE(q * mu)
        v, _ = Ld.eig(k=6)
        print(S, dt, 'K(E=q mu) eig:', np.round(v[np.argsort(-np.abs(v))], 4), flush=True)
