import sys, numpy as np
sys.path.insert(0, '.')
from ladder2 import Ladder
q = 4; mu = float(sys.argv[1]); dt = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
ch = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
Ld = Ladder(q, mu, 20, dt)
for E in np.round(np.linspace(0.02, 1.4, 24), 3):
    Ld.setE(E, channel=ch); v, _ = Ld.eig(k=6)
    v = v[np.argsort(-v.real)]
    print(f'E={E:.3f} E/mu={E/mu:6.2f}', ' '.join(f'{x.real:.4f}{x.imag:+.3f}j' for x in v[:5]), flush=True)
