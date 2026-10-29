import sys, numpy as np
sys.path.insert(0, '.')
from ladder_phys import Phys
from pole import pole
q = 4
mu = float(sys.argv[1])
z, _, _ = pole(q, mu); G1 = -z.real
print('2G1/mu', 2 * G1 / mu)
for S in [20, 40, 60]:
    ph = Phys(mu, S, 0.05)
    row = []
    for e in [3.0, 3.5, 3.8, 3.9, 3.95, 4.0, 4.05]:
        if e * mu < 2 * G1:
            row.append((e, round(ph.evals(e * mu)[0], 5)))
    print(S, row, flush=True)
