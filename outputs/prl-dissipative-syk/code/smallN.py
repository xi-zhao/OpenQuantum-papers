import sys, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian, dense_spectrum
for N in [6, 8]:
    for s in range(3):
        L = EvenLindbladian(N, seed=40 + s)
        row = []
        for mu in [0.005, 0.02, 0.1, 0.3]:
            ev = dense_spectrum(L, mu); r = np.sort(-ev.real)
            r = r[r > 1e-9]
            row.append(round(r[0] / mu, 3))
        print(N, s, 'gap/mu at mu=0.005,0.02,0.1,0.3:', row, flush=True)
