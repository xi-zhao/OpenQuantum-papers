"""Quasi-static heating test.  rho(0) = Gibbs(beta0).  Evolve rho with the Schrodinger
Lindbladian (even sector, D = -mu n).  Track E(t) (exact: E0 e^{-q mu t}) and the trace
distance between rho(t) and the Gibbs state rho_{beta(t)} with the same energy."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian
from scipy.optimize import brentq
N = int(sys.argv[1]); mu = float(sys.argv[2]); beta0 = float(sys.argv[3]); q = 4
L = EvenLindbladian(N, q=q, seed=int(sys.argv[4]) if len(sys.argv) > 4 else 0)
E, V = np.linalg.eigh(L.H)
def gibbs(b):
    w = np.exp(-b * (E - E.min())); w /= w.sum()
    return (V * w) @ V.conj().T
def energy(b):
    w = np.exp(-b * (E - E.min())); w /= w.sum(); return w @ E
def LS(v):
    A, B = L.split(v); nA, nB = L.size(A, B)
    return L.join(-1j * (L.Hp @ A - A @ L.Hp) - mu * nA, -1j * (L.Hm @ B - B @ L.Hm) - mu * nB)
def tnorm(v):
    A, B = L.split(v)
    return 0.5 * (np.abs(np.linalg.eigvalsh(A)).sum() + np.abs(np.linalg.eigvalsh(B)).sum())
rho = L.vec_of(gibbs(beta0)); h = L.vec_of(L.H)
E0 = (np.vdot(h, rho)).real
dt = 0.02; T = 3.0 / (q * mu); steps = int(T / dt); rec = []
for s in range(steps + 1):
    t = s * dt
    if s % max(1, steps // 60) == 0:
        Et = np.vdot(h, rho).real
        b = brentq(lambda x: energy(x) - Et, -1e-9, 200) if Et < -1e-12 else 0.0
        dG = tnorm(rho - L.vec_of(gibbs(b)))
        d0 = tnorm(L.vec_of(gibbs(beta0)) - L.vec_of(gibbs(b)))
        R = np.zeros((2*L.D,2*L.D),complex); R[np.ix_(L.ip,L.ip)], R[np.ix_(L.im,L.im)] = L.split(rho)
        Rd = V @ np.diag(np.diag(V.conj().T @ R @ V)) @ V.conj().T
        dD = tnorm(rho - L.vec_of(Rd))
        rec.append([t, Et / E0, np.exp(-q * mu * t), b, dG, d0, dD])
        print(*['%.5g' % x for x in rec[-1]], flush=True)
    k1 = LS(rho); k2 = LS(rho + 0.5 * dt * k1); k3 = LS(rho + 0.5 * dt * k2); k4 = LS(rho + dt * k3)
    rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
json.dump(rec, open(f'../data/therm_N{N}_mu{mu}_b{beta0}.json', 'w'))
