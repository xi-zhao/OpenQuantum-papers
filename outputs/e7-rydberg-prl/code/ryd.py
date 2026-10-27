"""Blockaded Rydberg chain: momentum-resolved exact diagonalization.

H = sum_j [ (1/2) X_j - Delta n_j ] + sum_{r>=2} V_r sum_j n_j n_{j+r},   Omega = 1,
constraint n_j n_{j+1} = 0 (nearest-neighbour blockade), V_r = V2 (2/r)^6.

Periodic chain of even length L; sectors k = 0 and k = pi.
Open chain (no translation) is also provided for finite-chain spectroscopy.
"""
import numpy as np
import numba as nb
import scipy.sparse as sp
import scipy.sparse.linalg as sla


# ---------------------------------------------------------------- basis
@nb.njit(cache=True)
def _count_bits(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@nb.njit(cache=True)
def _enum_constrained(L, periodic):
    # iterative DFS over bit strings with no two adjacent ones
    # number of states <= Fibonacci(L+2)
    a, b = 1, 2
    for _ in range(L - 1):
        a, b = b, a + b
    out = np.empty(b, dtype=np.int64)
    n = 0
    stack_state = np.empty(L + 2, dtype=np.int64)
    stack_pos = np.empty(L + 2, dtype=np.int64)
    # simple recursion emulation: build left to right
    # use explicit loop over integers instead (fast enough via Zeckendorf-style generation)
    top = 0
    stack_state[0] = 0
    stack_pos[0] = 0
    while top >= 0:
        s = stack_state[top]
        p = stack_pos[top]
        top -= 1
        if p == L:
            if periodic and L > 1 and (s & 1) and (s >> (L - 1)) & 1:
                continue
            out[n] = s
            n += 1
            continue
        # place 0
        top += 1
        stack_state[top] = s
        stack_pos[top] = p + 1
        # place 1 if previous bit empty
        if p == 0 or not ((s >> (p - 1)) & 1):
            top += 1
            stack_state[top] = s | (1 << p)
            stack_pos[top] = p + 1
    return np.sort(out[:n])


@nb.njit(cache=True)
def _rotl(s, L, mask):
    return ((s << 1) | (s >> (L - 1))) & mask


@nb.njit(cache=True)
def _rep_info(states, L):
    """For each state: representative (min over translations), shift s with state = T^s rep, period."""
    mask = (1 << L) - 1
    n = states.size
    rep = np.empty(n, dtype=np.int64)
    shift = np.empty(n, dtype=np.int64)
    for i in range(n):
        s0 = states[i]
        best = s0
        bs = 0
        t = s0
        for k in range(1, L):
            t = _rotl(t, L, mask)
            if t < best:
                best = t
                bs = k
        rep[i] = best
        # state = T^{-bs} applied to ... ; we store: rotl^bs(state) = rep  => state = rotr^bs(rep)
        shift[i] = bs
    return rep, shift


@nb.njit(cache=True)
def _period(r, L):
    mask = (1 << L) - 1
    t = _rotl(r, L, mask)
    p = 1
    while t != r:
        t = _rotl(t, L, mask)
        p += 1
    return p


class PeriodicBasis:
    def __init__(self, L):
        assert L % 2 == 0
        self.L = L
        st = _enum_constrained(L, True)
        rep, shift = _rep_info(st, L)
        self.reps = np.unique(rep)
        self.st_sorted = st
        self.st_rep_idx = np.searchsorted(self.reps, rep)
        self.st_shift = shift
        self.period = np.array([_period(r, L) for r in self.reps], dtype=np.int64)
        # occupation and staggered occupation of representatives
        bits = ((self.reps[:, None] >> np.arange(L)[None, :]) & 1).astype(np.int64)
        self.occ = bits.sum(1)
        self.stag = (bits * ((-1) ** np.arange(L))[None, :]).sum(1)
        # pair counts N_r = number of pairs at distance r (periodic, r = 2..L/2)
        self.pairs = np.zeros((self.reps.size, L // 2 + 1), dtype=np.int64)
        for d in range(2, L // 2 + 1):
            rolled = np.roll(bits, -d, axis=1)
            cnt = (bits * rolled).sum(1)
            if d == L // 2:
                cnt = cnt // 2  # each pair counted twice at antipodal distance
            self.pairs[:, d] = cnt


@nb.njit(cache=True)
def _build_offdiag(reps, period, st_sorted, st_rep_idx, st_shift, L, kpi):
    """Off-diagonal (X) part in momentum sector k = 0 (kpi=0) or pi (kpi=1).
    Basis |r,k> = 1/sqrt(R) sum_{s<R} e^{-iks} T^s |r>, valid if e^{ikR}=1."""
    n = reps.size
    rows = []
    cols = []
    vals = []
    valid = np.ones(n, dtype=np.bool_)
    if kpi:
        for i in range(n):
            if period[i] % 2 == 1:
                valid[i] = False
    for i in range(n):
        if not valid[i]:
            continue
        r = reps[i]
        for j in range(L):
            jm = (j - 1) % L
            jp = (j + 1) % L
            if (r >> j) & 1:
                c = r ^ (1 << j)
            else:
                if ((r >> jm) & 1) or ((r >> jp) & 1):
                    continue
                c = r | (1 << j)
            idx = np.searchsorted(st_sorted, c)
            ri = st_rep_idx[idx]
            if not valid[ri]:
                continue
            s = st_shift[idx]  # rotl^s(c) = rep  => c = T^{-s} rep (T = rotl)
            phase = 1.0
            if kpi and (s % 2 == 1):
                phase = -1.0
            v = 0.5 * phase * np.sqrt(period[i] / period[ri])
            rows.append(ri)
            cols.append(i)
            vals.append(v)
    return np.array(rows), np.array(cols), np.array(vals), valid


class PeriodicChain:
    def __init__(self, L):
        self.B = PeriodicBasis(L)
        self.L = L
        self._off = {}
        for kpi in (0, 1):
            B = self.B
            r, c, v, valid = _build_offdiag(B.reps, B.period, B.st_sorted, B.st_rep_idx,
                                            B.st_shift, L, kpi)
            keep = np.nonzero(valid)[0]
            pos = -np.ones(B.reps.size, dtype=np.int64)
            pos[keep] = np.arange(keep.size)
            X = sp.csr_matrix((v, (pos[r], pos[c])), shape=(keep.size, keep.size))
            self._off[kpi] = (X, keep)

    def H(self, Delta, V2, kpi, rmax=None):
        X, keep = self._off[kpi]
        B = self.B
        L = self.L
        if rmax is None:
            rmax = L // 2
        diag = -Delta * B.occ[keep].astype(float)
        for d in range(2, min(rmax, L // 2) + 1):
            diag = diag + V2 * (2.0 / d) ** 6 * B.pairs[keep, d]
        return X + sp.diags(diag)

    def ops(self, kpi):
        """Return vector representation of O_k acting from k=0 to sector kpi:
        O_0 = sum_j n_j (diagonal in k=0), O_pi = sum_j (-1)^j n_j maps k=0 -> k=pi."""
        B = self.B
        X0, keep0 = self._off[0]
        Xk, keepk = self._off[kpi]
        if kpi == 0:
            return B.occ[keep0].astype(float), None
        # O_pi |r,0> = stag_r |r,pi>  (only reps with even period survive; odd-period stag = 0)
        posk = -np.ones(B.reps.size, dtype=np.int64)
        posk[keepk] = np.arange(keepk.size)
        tgt = posk[keep0]
        coeff = B.stag[keep0].astype(float)
        return coeff, tgt


def lanczos_spectral(H, v0, m=200):
    """Lanczos with full reorthogonalization from v0. Returns Ritz energies and weights |<n|v0>|^2."""
    nrm = np.linalg.norm(v0)
    q = v0 / nrm
    n = v0.size
    Q = np.zeros((m, n))
    a = np.zeros(m)
    b = np.zeros(m)
    Q[0] = q
    w = H @ q
    a[0] = q @ w
    w -= a[0] * q
    for j in range(1, m):
        b[j] = np.linalg.norm(w)
        if b[j] < 1e-12:
            m = j
            break
        q = w / b[j]
        q -= Q[:j].T @ (Q[:j] @ q)
        q /= np.linalg.norm(q)
        Q[j] = q
        w = H @ q - b[j] * Q[j - 1]
        a[j] = q @ w
        w -= a[j] * q
    T = np.diag(a[:m]) + np.diag(b[1:m], 1) + np.diag(b[1:m], -1)
    e, U = np.linalg.eigh(T)
    wts = (nrm * U[0]) ** 2
    return e, wts


def ground(H):
    e, v = sla.eigsh(H, k=1, which='SA', tol=1e-12)
    return e[0], v[:, 0]


def lowest(H, k=6):
    e = sla.eigsh(H, k=k, which='SA', tol=1e-10, return_eigenvectors=False)
    return np.sort(e)


# ---------------------------------------------------------------- open chain
class OpenChain:
    def __init__(self, L):
        self.L = L
        st = _enum_constrained(L, False)
        self.st = st
        bits = ((st[:, None] >> np.arange(L)[None, :]) & 1).astype(np.int64)
        self.bits = bits
        rows, cols = [], []
        for j in range(L):
            flip = st ^ (1 << j)
            ok = np.ones(st.size, bool)
            if j > 0:
                ok &= ((flip >> (j - 1)) & 1) == 0
            if j < L - 1:
                ok &= ((flip >> (j + 1)) & 1) == 0
            idx = np.searchsorted(st, flip[ok])
            rows.append(idx)
            cols.append(np.nonzero(ok)[0])
        rows = np.concatenate(rows)
        cols = np.concatenate(cols)
        self.X = sp.csr_matrix((0.5 * np.ones(rows.size), (rows, cols)), shape=(st.size, st.size))
        self.pairs = np.zeros((st.size, L), dtype=np.int64)
        for d in range(2, L):
            self.pairs[:, d] = (bits[:, :-d] * bits[:, d:]).sum(1)

    def H(self, Delta, V2, local=None):
        diag = -Delta * self.bits.sum(1).astype(float)
        for d in range(2, self.L):
            diag = diag + V2 * (2.0 / d) ** 6 * self.pairs[:, d]
        if local is not None:
            diag = diag - self.bits @ local
        return self.X + sp.diags(diag)
