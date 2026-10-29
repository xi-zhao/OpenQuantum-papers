"""Check SD at mu>0: Heisenberg evolve psi_i (odd sector) at finite N, compare C(t)=Tr[psi_i(t) psi_i]/d."""
import sys, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian
from sd import solve
N = int(sys.argv[1]); mu = float(sys.argv[2])
L = EvenLindbladian(N, seed=11)
D = L.D; Hp, Hm = L.Hp, L.Hm; c = [ci.toarray() for ci in L.c]
def apply(C, Dm):
    SC = sum(ci @ Dm @ ci for ci in c); SD = sum(ci.conj().T @ C @ ci.conj().T for ci in c)
    return (1j * (Hp @ C - C @ Hm) - mu * (0.5 * N * C + SC),
            1j * (Hm @ Dm - Dm @ Hp) - mu * (0.5 * N * Dm + SD))
acc = None
for i in range(3):
    C0, D0 = c[i].copy(), c[i].conj().T.copy()
    C, Dm = C0.copy(), D0.copy(); dt = 0.02; out = []
    for s in range(301):
        if s % 25 == 0:
            out.append((np.vdot(C0, C) + np.vdot(D0, Dm)).real / (2 * D))
        k1 = apply(C, Dm); k2 = apply(C + dt/2*k1[0], Dm + dt/2*k1[1])
        k3 = apply(C + dt/2*k2[0], Dm + dt/2*k2[1]); k4 = apply(C + dt*k3[0], Dm + dt*k3[1])
        C = C + dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); Dm = Dm + dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    acc = np.array(out) if acc is None else acc + np.array(out)
acc /= 3
t, Csd = solve(q=4, mu=mu, T=6, dt=0.002)
for k, v in enumerate(acc):
    tt = k * 0.5; print(tt, round(v, 4), round(np.interp(tt, t, Csd), 4))
