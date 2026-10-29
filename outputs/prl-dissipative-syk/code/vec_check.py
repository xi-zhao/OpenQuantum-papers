"""Check vectorized two-site Green functions at finite N:
G_LL(t>0)=C, G_LR(t>0)=iC with R_i(O)=i P O P psi_i, L_i(O)=psi_i O."""
import sys, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian
N = 8
L = EvenLindbladian(N, seed=3)
d = 2 ** (N // 2)
H = L.H; psi = [p.toarray() for p in L.psi]
par = np.diag([(-1) ** bin(b).count('1') for b in range(d)]).astype(complex)
mu = 0.3
def Lsup(O):
    n = 0.5 * N * O + 1j * sum(p @ (1j * par @ O @ par @ p) for p in psi)  # N/2 + i sum L R
    return 1j * (H @ O - O @ H) - mu * n
# check n on psi_0 (odd) -> size 1, and on psi0 psi1 -> 2
for O, s in [(psi[0], 1), (psi[0] @ psi[1], 2), (psi[0] @ psi[1] @ psi[2] @ psi[3], 4)]:
    nO = 0.5 * N * O + 1j * sum(p @ (1j * par @ O @ par @ p) for p in psi)
    print('size', s, np.allclose(nO, s * O))
