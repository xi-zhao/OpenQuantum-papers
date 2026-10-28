"""Float64 PIQS local-emission generator, in the existing Dicke-block convention.

The gamma1--gamma4 formulas are those of QuTiP PIQS (_piqs.pyx,
QuTiP 5.3.1; Shammah et al., PRA 98, 063815), evaluated in float64.
Only independent local emission with rate one is supported here. This does
not patch QuTiP or change the generator used by earlier research scripts.
"""
from functools import lru_cache
import numpy as np
import scipy.sparse as sp
from qutip import piqs


@lru_cache(maxsize=2)
def space(N):
  jx, _, jz = piqs.jspin(N, basis='dicke')
  ne = np.rint(jz.diag().real + N / 2).astype(int)
  twoj = np.concatenate([np.full(j + 1, j, dtype=int)
                         for j in range(N, -1, -2)])
  twom = 2 * ne - N
  mask = piqs.block_matrix(N).toarray().astype(bool)
  idx = np.flatnonzero(mask.ravel(order='F'))
  d = len(ne)
  row, col = idx % d, idx // d
  lookup = {(int(twoj[i]), int(twom[i]), int(twom[j])): k
            for k, (i, j) in enumerate(zip(row, col))}
  rows, cols, vals = [], [], []
  for r, (i, j) in enumerate(zip(row, col)):
    tj, tm, tm1 = int(twoj[i]), int(twom[i]), int(twom[j])
    rows.append(r); cols.append(r); vals.append(-(N + (tm + tm1)/2)/2)
    for dj in (0, 2, -2):
      src = (tj + dj, tm + 2, tm1 + 2)
      c = lookup.get(src)
      if c is None:
        continue
      a, m, m1 = (x / 2 for x in src)
      if dj == 0:
        value = (0.5 * np.sqrt((a+m)*(a-m+1)*(a+m1)*(a-m1+1))
                 * (N/2+1)/(a*(a+1))) if a > 0 else 0.0
      elif dj == 2:
        value = (0.5 * np.sqrt((a+m)*(a+m-1)*(a+m1)*(a+m1-1))
                 * (N/2+a+1)/(a*(2*a+1))) if a > 0 else 0.0
      else:
        value = (0.5 * np.sqrt((a-m+1)*(a-m+2)*(a-m1+1)*(a-m1+2))
                 * (N/2-a)/((a+1)*(2*a+1)))
      rows.append(r); cols.append(c); vals.append(value)
  decay = sp.csc_matrix((vals, (rows, cols)), shape=(len(idx), len(idx)))
  swap = np.empty(len(idx), dtype=int)
  for k, (i, j) in enumerate(zip(row, col)):
    swap[k] = lookup[(int(twoj[j]), int(twom[j]), int(twom[i]))]
  return dict(N=N, d=d, idx=idx, row=row, col=col, ne=ne,
              twoj=twoj, jx=sp.csc_matrix(jx.data.as_scipy()),
              decay=decay, swap=swap, trace=(row == col).astype(float))


def generators(N, omega=2.0, interaction=10.0, backend='float64'):
  s = space(N)
  ne, d, idx = s['ne'], s['d'], s['idx']
  h0 = omega * s['jx'] + sp.diags(interaction/(2*N)*ne*(ne-1))
  eye = sp.eye(d, format='csc')
  comm = -1j * (sp.kron(eye, h0, format='csc')
                - sp.kron(h0.T, eye, format='csc'))
  if backend == 'piqs':
    raw = piqs.Dicke(N=N, emission=1.0).lindbladian().data.as_scipy()
    decay = sp.csc_matrix(raw)[idx][:, idx]
  elif backend == 'float64':
    decay = s['decay']
  else:
    raise ValueError(backend)
  L0 = (comm[idx][:, idx] + decay).tocsc()
  derivative = sp.diags(1j*(ne[s['row']] - ne[s['col']]), format='csc')
  return L0, derivative, s


def populations(vec, s):
  mask = s['row'] == s['col']
  k = s['ne'][s['row'][mask]]
  return (np.bincount(k, weights=vec[mask].real, minlength=s['N']+1)
          + 1j*np.bincount(k, weights=vec[mask].imag, minlength=s['N']+1))
