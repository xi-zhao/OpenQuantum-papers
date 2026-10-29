r"""Second-moment closure check: B = sum_i [psi_i,H]^dag [psi_i,H].
L(H^2) = -2 q mu H^2 + mu B.  If B -> b*1 (relative fluctuation -> 0), H^2 - c is an
approximate eigenoperator with rate 2 q mu.  Report std(B)/mean(B) and the exact rate of
the k=2 secular mode."""
import sys, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian
for N in [int(x) for x in sys.argv[1].split(',')]:
    rs = []
    for s in range(3):
        L = EvenLindbladian(N, seed=9000 + s)
        H = L.H; d = H.shape[0]
        B = np.zeros_like(H)
        for p in L.psi:
            p = p.toarray(); X = p @ H - H @ p; B += X.conj().T @ X
        m = np.trace(B).real / d
        dev = B - m * np.eye(d)
        rs.append([m, np.sqrt(np.trace(dev.conj().T @ dev).real / d) / m])
    rs = np.array(rs).mean(0)
    print(N, 'mean B', round(rs[0], 4), ' rel.fluct', round(rs[1], 4), ' 2q sigma^2', round(2 * 4 * np.trace(H @ H).real / d, 4), flush=True)
