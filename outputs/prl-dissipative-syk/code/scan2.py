import sys, time, json, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian, leading
N = int(sys.argv[1]); seeds = int(sys.argv[2])
mus = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8]
out = []
for s in range(seeds):
    L = EvenLindbladian(N, seed=100 + s)
    v0 = None
    for mu in mus:
        t = time.time()
        v, V = leading(L, mu, k=3, tol=1e-7)
        rec = dict(N=N, seed=s, mu=mu, lam=[[x.real, x.imag] for x in v],
                   size=float(L.mean_size(V[:, 0])))
        out.append(rec)
        print(N, s, mu, f'{-v[0].real:.4f} {-v[0].real/mu:.3f} im={v[0].imag:.3f}', f'{time.time()-t:.1f}s', flush=True)
json.dump(out, open(f'../data/gap_N{N}.json', 'w'))
