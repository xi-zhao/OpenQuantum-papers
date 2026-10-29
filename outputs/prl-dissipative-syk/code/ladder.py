r"""Large-N bound states of the vectorized Lindbladian (two-site 'Euclidean' T=0 picture).
-Lc = Hhat = -i H_L + i H_R + mu N/2 + i mu sum psi_L psi_R,  ground state |1>>.
G_LL=G_RR=sgn(t)C(|t|), G_LR=iC(|t|), G_RL=-iC(|t|); C from Kadanoff-Baym (sd.py).
Kernel (MS convention): K_{ab,cd} = -(q-1) J^2 c_c c_d G_ac(t13) G_bd(t24) G_cd(t34)^{q-2},
c_L=-i, c_R=+i.  Bound state of decay rate E: Psi = e^{-E(t1+t2)/2} phi(t12), K phi = phi."""
import sys, numpy as np
sys.path.insert(0, '.')
from sd import solve
from scipy.signal import fftconvolve
from scipy.sparse.linalg import eigs

FL = ['L', 'R']
cfac = {'L': -1j, 'R': 1j}

def greens(q, mu, S, dt, J=1.0):
    t, C = solve(q=q, J=J, mu=mu, T=S + dt, dt=dt)
    n = int(round(S / dt))
    s = np.arange(-n, n + 1) * dt
    Ca = np.interp(np.abs(s), t, C)
    sg = np.sign(s)
    G = {('L', 'L'): sg * Ca, ('R', 'R'): sg * Ca,
         ('L', 'R'): 1j * Ca, ('R', 'L'): -1j * Ca}
    return s, G

def kernel(q, mu, E, S=30.0, dt=0.05, J=1.0, channel=1.0, G=None):
    """channel: 1 for singlet; +-1/(q-1) for non-singlet."""
    if G is None:
        s, G = greens(q, mu, S, dt, J)
    else:
        s, G = G
    n = len(s)
    ew = np.exp(E * s / 2)
    g = {k: v * ew for k, v in G.items()}
    blocks = np.zeros((4 * n, 4 * n), complex)
    idx = {(a, b): i for i, (a, b) in enumerate([(a, b) for a in FL for b in FL])}
    # M_{ab,cd}(x) = int du g_ac(u) g_bd(u - x)  -> correlation; x on grid of size 2n-1
    for (a, b), I in idx.items():
        for (c, d), Jx in idx.items():
            Mx = fftconvolve(g[(a, c)], g[(b, d)][::-1], mode='full') * dt  # index k -> x=(k-(n-1))*dt
            w = -(q - 1) * J**2 * cfac[c] * cfac[d] * G[(c, d)] ** (q - 2) * channel
            # K(s_i, r_j) = M(s_i - r_j) w(r_j) dt
            ii = np.arange(n)[:, None]; jj = np.arange(n)[None, :]
            blocks[I * n:(I + 1) * n, Jx * n:(Jx + 1) * n] = Mx[ii - jj + n - 1] * w[None, :] * dt
    return s, blocks

def near_one(K, k=4):
    vals = eigs(K, k=k, sigma=1.0, return_eigenvectors=False)
    return vals[np.argsort(np.abs(vals - 1))]

if __name__ == '__main__':
    q = 4; mu = float(sys.argv[1])
    S = float(sys.argv[2]) if len(sys.argv) > 2 else 30.0
    dt = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
    Gs = greens(q, mu, S, dt)
    for E in np.linspace(0.2 * q * mu, min(2.4 * q * mu, 1.2), 12):
        s, K = kernel(q, mu, E, G=Gs)
        v = near_one(K)
        print(f'E={E:.4f} E/(q mu)={E/(q*mu):.3f}', np.round(v[:3], 4), flush=True)
