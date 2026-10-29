"""Odd-operator sector (governs fermion Green's functions).  Keldysh-consistent
dissipator -mu*n with n(string)=size; on odd O: n(O) = N/2 O + sum_i psi_i O psi_i.
Odd O = [[0, C], [D, 0]] in (+,-) parity blocks."""
import sys, json, numpy as np
sys.path.insert(0, '.')
from syk_ed import EvenLindbladian
from scipy.sparse.linalg import LinearOperator, eigs
N = int(sys.argv[1]); seed = int(sys.argv[2])
L = EvenLindbladian(N, seed=7000 + seed)
D = L.D; Hp, Hm = L.Hp, L.Hm
c = [ci.toarray() for ci in L.c]
def apply(v, mu):
    C = v[:D*D].reshape(D, D); Dm = v[D*D:].reshape(D, D)
    # psi = [[0,c],[c^dag,0]];  psi O psi for O=[[0,C],[Dm,0]] = [[0, c Dm c], [c^dag C c^dag, 0]]
    SC = sum(ci @ Dm @ ci for ci in c); SD = sum(ci.conj().T @ C @ ci.conj().T for ci in c)
    nC = 0.5 * N * C + SC; nD = 0.5 * N * Dm + SD
    LC = 1j * (Hp @ C - C @ Hm) - mu * nC
    LD = 1j * (Hm @ Dm - Dm @ Hp) - mu * nD
    return np.concatenate([LC.ravel(), LD.ravel()])
# sanity: size of psi_0 must be 1
p0 = L.psi[0].toarray(); v = np.concatenate([p0[np.ix_(L.ip, L.im)].ravel(), p0[np.ix_(L.im, L.ip)].ravel()])
r = apply(v, 1.0); print('size(psi0)=', -np.vdot(v, r).real / np.vdot(v, v).real)
out = []
for mu in [0.0, 0.1, 0.2, 0.4, 0.8]:
    op = LinearOperator((2*D*D,)*2, matvec=lambda x: apply(x, mu), dtype=complex)
    w = eigs(op, k=4, which='LR', tol=1e-7, ncv=40, maxiter=20000, return_eigenvectors=False)
    w = w[np.argsort(-w.real)]
    out.append([mu, -w[0].real, abs(w[0].imag)]); print(N, seed, out[-1], flush=True)
json.dump(out, open(f'../data/odd_N{N}_s{seed}.json', 'w'))
