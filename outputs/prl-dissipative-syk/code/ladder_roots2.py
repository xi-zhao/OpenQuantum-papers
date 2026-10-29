import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder2 import Ladder
from pole import pole
from scipy.optimize import brentq
q = 4
def evs(Ld, E, k=8, ch=1.0):
    Ld.setE(E, channel=ch); v, _ = Ld.eig(k=k)
    return v
def track(mu, S=20, dt=0.05, nE=30):
    z, _, _ = pole(q, mu); G1 = -z.real
    Ld = Ladder(q, mu, S, dt)
    Es = np.linspace(0.005, 0.99 * 2 * G1, nE)
    data = []
    for E in Es:
        v = evs(Ld, E); r = np.sort(v[np.abs(v.imag) < 1e-6].real)[::-1]
        data.append(r)
    roots = []
    def g(E):
        v = evs(Ld, E); r = v[np.abs(v.imag) < 1e-6].real
        return r[np.argmin(np.abs(r - 1))] - 1
    for i in range(nE - 1):
        a = np.sum(data[i] > 1); b = np.sum(data[i + 1] > 1)
        if a != b:
            try: roots.append(brentq(g, Es[i], Es[i + 1], xtol=1e-6))
            except ValueError: roots.append(np.nan)
    return G1, roots, Es, data
if __name__ == '__main__':
    out = []
    for mu in [float(x) for x in sys.argv[1].split(',')]:
        G1, roots, Es, data = track(mu)
        print(f'mu={mu} 2G1/mu={2*G1/mu:.3f} roots E/mu=', [round(r / mu, 4) for r in roots],
              ' n(>1) at E->0:', np.sum(data[0] > 1), flush=True)
        out.append(dict(mu=mu, G1=G1, roots=roots))
    json.dump(out, open('../data/ladder_roots_%s.json' % sys.argv[2], 'w'))
