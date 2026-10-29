import sys, numpy as np
sys.path.insert(0, '.')
from ladder2 import Ladder
q = 4; mu = float(sys.argv[1]); ch = float(sys.argv[2]); dt = float(sys.argv[3]); S = float(sys.argv[4])
Ld = Ladder(q, mu, S, dt)
for E in np.round(np.arange(0.0, 2.2 * mu + 1e-9, 0.2 * mu), 5):
    Ld.setE(E, channel=ch); v, _ = Ld.eig(k=8)
    r = np.sort(v.real[np.abs(v.imag) < 1e-6])[::-1]
    print(f'E/mu={E/mu:4.1f}', np.round(r[:6], 4))
