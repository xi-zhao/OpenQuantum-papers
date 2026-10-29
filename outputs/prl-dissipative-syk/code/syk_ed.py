"""Exact diagonalization of the dissipative-SYK Lindbladian on parity-even operators.

Majoranas: {psi_i, psi_j} = delta_ij (psi^2 = 1/2).  Jump operators L_i = sqrt(mu) psi_i.
On even operators the adjoint (Heisenberg) Lindbladian is exactly
    L(O) = i[H, O] - mu * n(O),   n(O) = (N/2) O - sum_i psi_i O psi_i,
where n multiplies a Majorana string of length s by s (operator size).
Even operators are block diagonal O = diag(A, B) in the (+,-) parity basis.
"""
import itertools
from math import comb

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, eigs


def majoranas(N):
    """Jordan-Wigner Majoranas, normalized so psi^2 = 1/2."""
    n = N // 2
    X = sp.csr_matrix([[0, 1], [1, 0]], dtype=complex)
    Y = sp.csr_matrix([[0, -1j], [1j, 0]])
    Z = sp.csr_matrix([[1, 0], [0, -1]], dtype=complex)
    I = sp.identity(2, dtype=complex, format="csr")
    out = []
    for k in range(n):
        for P in (X, Y):
            M = None
            for o in [Z] * k + [P] + [I] * (n - k - 1):
                M = o if M is None else sp.kron(M, o, format="csr")
            out.append(M / np.sqrt(2))
    return out


class EvenLindbladian:
    def __init__(self, N, q=4, J=1.0, seed=0):
        assert q % 4 == 0, "Hermitian without phase only for q = 0 mod 4"
        self.N, self.q, self.J = N, q, J
        rng = np.random.default_rng(seed)
        psi = majoranas(N)
        d = 2 ** (N // 2)
        par = np.array([bin(b).count("1") % 2 for b in range(d)])
        self.ip, self.im = np.where(par == 0)[0], np.where(par == 1)[0]
        # exact normalization: <[H,psi_i]^2> = 2^{1-q} J^2 at every N
        sigma = np.sqrt(J**2 / comb(N - 1, q - 1))
        H = sp.csr_matrix((d, d), dtype=complex)
        for idx in itertools.combinations(range(N), q):
            M = psi[idx[0]]
            for i in idx[1:]:
                M = M @ psi[i]
            H = H + rng.normal(0.0, sigma) * M
        H = H.toarray()
        self.H = 0.5 * (H + H.conj().T)
        self.Hp = self.H[np.ix_(self.ip, self.ip)]
        self.Hm = self.H[np.ix_(self.im, self.im)]
        # psi_i = [[0, c_i], [c_i^dag, 0]] in the (+,-) basis
        self.c = [p[self.ip][:, self.im].tocsr() for p in psi]
        self.cdT = [ci.conj().tocsr() for ci in self.c]  # (c^dag)^T = conj(c)
        self.D = d // 2
        self.psi = psi

    def split(self, v):
        D = self.D
        return v[: D * D].reshape(D, D), v[D * D :].reshape(D, D)

    def join(self, A, B):
        return np.concatenate([A.ravel(), B.ravel()])

    def size(self, A, B):
        SA = np.zeros_like(A)
        SB = np.zeros_like(B)
        for c, cT in zip(self.c, self.cdT):
            # c B c^dag  and  c^dag A c
            SA += (cT @ (c @ B).T).T
            SB += (c.T @ (c.conj().T @ A).T).T
        return 0.5 * self.N * A - SA, 0.5 * self.N * B - SB

    def apply(self, v, mu):
        A, B = self.split(v)
        nA, nB = self.size(A, B)
        LA = 1j * (self.Hp @ A - A @ self.Hp) - mu * nA
        LB = 1j * (self.Hm @ B - B @ self.Hm) - mu * nB
        return self.join(LA, LB)

    def vec_of(self, O):
        return self.join(O[np.ix_(self.ip, self.ip)], O[np.ix_(self.im, self.im)])

    def mean_size(self, v):
        A, B = self.split(v)
        nA, nB = self.size(A, B)
        return (np.vdot(v, self.join(nA, nB)) / np.vdot(v, v)).real


def deflation_basis(L, include_H=True):
    one = L.vec_of(np.eye(2 * L.D))
    one = one / np.linalg.norm(one)
    basis = [one]
    if include_H:
        h = L.vec_of(L.H)
        h = h - np.vdot(one, h) * one
        basis.append(h / np.linalg.norm(h))
    return np.array(basis)


def leading(L, mu, k=6, include_H=True, tol=1e-8, v0=None):
    """Eigenvalues with largest real part on the complement of span{1, H}.

    1 and H are simultaneously left and right eigenvectors of L, so the
    orthogonal complement is invariant and the deflation is exact."""
    n = 2 * L.D**2
    B = deflation_basis(L, include_H)
    shift = 10.0 * mu * L.N + 50.0

    def P(v):
        return v - B.T @ (B.conj() @ v)

    def mv(v):
        w = P(v)
        return P(L.apply(w, mu)) - shift * (v - w)

    op = LinearOperator((n, n), matvec=mv, dtype=complex)
    vals, vecs = eigs(op, k=k, which="LR", tol=tol, ncv=max(4 * k, 40), v0=v0,
                      maxiter=20000)
    o = np.argsort(-vals.real)
    return vals[o], vecs[:, o]


def dense_spectrum(L, mu):
    n = 2 * L.D**2
    M = np.empty((n, n), dtype=complex)
    e = np.zeros(n, dtype=complex)
    for j in range(n):
        e[j] = 1.0
        M[:, j] = L.apply(e, mu)
        e[j] = 0.0
    return np.linalg.eigvals(M)
