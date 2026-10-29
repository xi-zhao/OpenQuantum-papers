import sys, json, time, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian, leading
N = int(sys.argv[1]); seed = int(sys.argv[2])
L = EvenLindbladian(N, seed=2000 + seed)
out = []
for mu in [float(x) for x in (sys.argv[3] if len(sys.argv)>3 else "0.8,0.6,0.4,0.3,0.2").split(",")]:
    t = time.time()
    v, V = leading(L, mu, k=2, tol=1e-6)
    out.append([mu, -v[0].real, v[0].imag, L.mean_size(V[:, 0])])
    print(N, seed, mu, out[-1], '%.0fs' % (time.time() - t), flush=True)
    json.dump(out, open(f'../data/fixedmu_N{N}_s{seed}.json', 'w'))
