import sys, numpy as np
sys.path.insert(0, '.')
from ladder2 import Ladder, PAIRS
q = 4; mu = float(sys.argv[1])
Ld = Ladder(q, mu, 20, 0.05); n = Ld.n
def X(v):
    out = np.zeros_like(v)
    for i, (a, b) in enumerate(PAIRS):
        j = PAIRS.index((b, a))
        out[i*n:(i+1)*n] = -v[j*n:(j+1)*n][::-1]
    return out
for E in [float(x) * mu for x in sys.argv[2].split(',')]:
    Ld.setE(E); v, V = Ld.eig(k=8)
    o = np.argsort(-v.real)
    for k in o[:6]:
        x = V[:, k]; par = np.vdot(x, X(x)).real / np.vdot(x, x).real
        print(f'E/mu={E/mu:.2f} lam={v[k].real:.4f}{v[k].imag:+.1e}j  X-parity={par:+.3f}')
