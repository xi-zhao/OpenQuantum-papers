import sys, json, time, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian, leading
N = int(sys.argv[1]); seeds = [int(s) for s in sys.argv[2].split(',')]
mus = [0.8, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05]
out = []
for s in seeds:
    L = EvenLindbladian(N, seed=5000 + s)
    for mu in mus:
        v, V = leading(L, mu, k=2, tol=1e-6)
        out.append([s, mu, -v[0].real, v[0].imag, L.mean_size(V[:, 0])])
        print(*out[-1], flush=True)
        json.dump(out, open(f'../data/grid_N{N}.json', 'w'))
