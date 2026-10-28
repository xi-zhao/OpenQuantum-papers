"""Small-N validation against a full local-spin Lindblad construction.

The reference is projected onto orthonormal symmetrized Pauli products, an
independent representation of the entire permutation-invariant operator space.
"""
from itertools import permutations, product
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigvals
from scipy.optimize import linear_sum_assignment
from biorthogonal_generator import generators


def mm(a,b):
  return np.einsum('ik,kj->ij',a,b,optimize=False)


def kron_all(items):
  out=np.ones((1,1),complex)
  for a in items:
    out=np.kron(out,a)
  return out


def validate(N):
  one=np.eye(2,dtype=complex)
  sx=np.array([[0,1],[1,0]],complex)
  sy=np.array([[0,-1j],[1j,0]],complex)
  sz=np.diag([1.,-1.])
  n=np.diag([1.,0.])
  jump=np.array([[0,0],[1,0]],complex)
  local=lambda op,i:kron_all([op if j==i else one for j in range(N)])
  Ne=sum(local(n,i) for i in range(N))
  H=-3.25*Ne+sum(local(sx,i) for i in range(N))+5/N*(mm(Ne,Ne)-Ne)
  eye=np.eye(2**N)
  L=-1j*(np.kron(eye,H)-np.kron(H.T,eye))
  for i in range(N):
    J=local(jump,i); JJ=mm(J.conj().T,J)
    L+=np.kron(J.conj(),J)-0.5*(np.kron(eye,JJ)+np.kron(JJ.T,eye))
  pauli=[one,sx,sy,sz]
  columns=[]
  for counts in product(range(N+1),repeat=3):
    if sum(counts)>N:
      continue
    counts=(*counts,N-sum(counts))
    order=tuple(a for a in range(4) for _ in range(counts[a]))
    terms=set(permutations(order))
    mat=sum(kron_all([pauli[a] for a in term]) for term in terms)
    columns.append(mat.ravel(order='F')/np.sqrt(len(terms)*2**N))
  Q=np.column_stack(columns)
  A=np.einsum('ki,kj->ij',Q.conj(),mm(L,Q),optimize=False)
  L0,L1,s=generators(N)
  B=(L0+3.25*L1).toarray()
  wa,wb=eigvals(A),eigvals(B)
  i,j=linear_sum_assignment(abs(wa[:,None]-wb[None,:]))
  out=dict(N=N,dimension=len(wb),spectral_max_error=float(max(abs(wa[i]-wb[j]))),
           reference_invariance_error=float(np.linalg.norm(mm(L,Q)-mm(Q,A))),
           trace_defect=float(max(abs(np.einsum('i,ij->j',s['trace'],B,optimize=False)))),
           hermiticity_defect=float(np.max(abs(B-B[s['swap']][:,s['swap']].conj()))))
  assert out['spectral_max_error']<1e-10, out
  assert out['reference_invariance_error']<1e-11, out
  assert out['trace_defect']<1e-12, out
  assert out['hermiticity_defect']<1e-12, out
  return out


if __name__=='__main__':
  result=[validate(N) for N in (1,2,3,4)]
  path=Path(__file__).resolve().parents[1]/'data/biorthogonal_audit'
  path.mkdir(exist_ok=True)
  (path/'generator_validation.json').write_text(json.dumps(result,indent=2)+'\n')
  print(json.dumps(result,indent=2))
