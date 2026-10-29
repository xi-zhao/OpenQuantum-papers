import sys, time, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian, leading
N = int(sys.argv[1]); seeds = int(sys.argv[2])
mus = [0.02, 0.05, 0.1, 0.2, 0.4, 0.8]
for s in range(seeds):
    L = EvenLindbladian(N, seed=s)
    for mu in mus:
        t = time.time()
        v, V = leading(L, mu, k=4)
        print(N, s, mu, ' '.join(f'{x.real/mu:.3f}{x.imag/mu:+.3f}j' for x in v[:4]),
              f'size={L.mean_size(V[:,0]):.2f}', f'{time.time()-t:.1f}s', flush=True)
