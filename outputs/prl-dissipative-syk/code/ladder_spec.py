r"""Bound states of the large-N ladder: real E where K(E) has eigenvalue 1 (singlet)
or -(q-1) (non-singlet, since K_ns = -K/(q-1)).  Dense kernel from ladder.py."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from ladder import greens, kernel
from pole import pole
from scipy.linalg import eigvals
from scipy.optimize import brentq
q = 4

def kev(mu, E, Gs):
    s, K = kernel(q, mu, E, G=Gs)
    return eigvals(K)

def count(v, thr, side):
    r = v[np.abs(v.imag) < 1e-7].real
    return int(np.sum(r > thr)) if side > 0 else int(np.sum(r < thr))

def closest(mu, E, Gs, thr):
    v = kev(mu, E, Gs); r = v[np.abs(v.imag) < 1e-7].real
    return r[np.argmin(np.abs(r - thr))] - thr

def bound_states(mu, S=20.0, dt=0.1, nE=40):
    z, _, _ = pole(q, mu); G1 = -z.real
    Emax = 2 * G1 * 0.985
    Gs = greens(q, mu, S, dt)
    Es = np.linspace(0.01, Emax, nE)
    cs = {'singlet': [], 'nonsinglet': []}
    for E in Es:
        v = kev(mu, E, Gs)
        cs['singlet'].append(count(v, 1.0, +1)); cs['nonsinglet'].append(count(v, -(q - 1), -1))
    out = {'mu': mu, 'Gamma1': G1, 'singlet': [], 'nonsinglet': []}
    for ch, thr in [('singlet', 1.0), ('nonsinglet', -(q - 1.0))]:
        c = cs[ch]
        for i in range(nE - 1):
            if c[i] != c[i + 1]:
                try:
                    r = brentq(lambda E: closest(mu, E, Gs, thr), Es[i], Es[i + 1], xtol=1e-5)
                except ValueError:
                    r = 0.5 * (Es[i] + Es[i + 1])
                out[ch].append(r)
    return out

if __name__ == '__main__':
    res = []
    for mu in [float(x) for x in sys.argv[1].split(',')]:
        o = bound_states(mu)
        res.append(o)
        print(f"mu={mu} 2G1={2*o['Gamma1']:.4f} singlet E/mu={[round(e/mu,3) for e in o['singlet']]} nonsinglet E/mu={[round(e/mu,3) for e in o['nonsinglet']]}", flush=True)
    json.dump(res, open(f'../data/ladder_{sys.argv[2] if len(sys.argv)>2 else "scan"}.json', 'w'))
