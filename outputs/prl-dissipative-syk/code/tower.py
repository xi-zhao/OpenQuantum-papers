"""Energy tower: orthogonal polynomials P_k(H) (w.r.t. infinite-T trace) and their
Majorana size.  Prediction: n(P_k) ~ k q P_k (OU / q-OU tower), up to 1/N corrections."""
import sys, numpy as np, itertools
sys.path.insert(0, '.')
from syk_ed import majoranas
from math import comb
N = int(sys.argv[1]); q = 4; K = 5
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
psi = [p.toarray() for p in majoranas(N)]
d = psi[0].shape[0]
H = np.zeros((d, d), complex)
sig = np.sqrt(1 / comb(N - 1, q - 1))
for idx in itertools.combinations(range(N), q):
    M = psi[idx[0]]
    for i in idx[1:]:
        M = M @ psi[i]
    H += rng.normal(0, sig) * M
H = (H + H.conj().T) / 2
E, V = np.linalg.eigh(H)
ip = lambda A, B: np.vdot(A, B).real / d
# orthogonal polynomials via Lanczos in energy eigenbasis (diagonal)
P = [np.ones(d)]
x = E.copy()
for k in range(1, K + 1):
    p = x * P[-1]
    for r in P:
        p = p - (r @ p) / (r @ r) * r
    P.append(p)
psiE = [V.conj().T @ p @ V for p in psi]
def size(Ddiag):
    O = np.diag(Ddiag).astype(complex)
    S = 0.5 * N * O
    for pe in psiE:
        S = S - pe @ O @ pe
    return S
# eta (chord intersection weight) for Majorana SYK
eta = sum((-1) ** k * comb(q, k) * comb(N - q, q - k) for k in range(q + 1)) / comb(N, q)
print('N', N, 'eta', round(eta, 4))
for k in range(1, K + 1):
    S = size(P[k])
    Pk = np.diag(P[k])
    s = np.vdot(Pk, S).real / np.vdot(Pk, Pk).real
    diagS = np.diag(np.diag(S))
    # component that stays a function of H
    leak = np.linalg.norm(S - diagS) / np.linalg.norm(S)
    # component within functions of H but different polynomial
    rest = np.diag(S).real - s * P[k]
    mix = np.linalg.norm(rest) / np.linalg.norm(S)
    print(f'k={k} size/q={s/q:.4f}  offdiag(leak)={leak:.4f}  diag-mix={mix:.4f}', flush=True)
