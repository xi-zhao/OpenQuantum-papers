# full hybrid Liouvillian, n sites PBC ring, t=1,gamma=1, loss L_j=c_j (post-selected: -n gamma rho), perturb z*sum L_xi rho L_xi^dag, L_xi = c^dag_{xi1}..c^dag_{xid}
import numpy as np, itertools, scipy.sparse as sp, scipy.sparse.linalg as sla, sys
from pblock import Nmat
def ops(n):
    I=sp.identity(2,format='csr'); Z=sp.csr_matrix(np.diag([1.,-1.])); a=sp.csr_matrix(np.array([[0.,1.],[0.,0.]]))  # |1>->|0>, basis (0,1)? a|1>=|0>
    cs=[]
    for j in range(n):
        m=[Z]*j+[a]+[I]*(n-j-1); o=m[0]
        for x in m[1:]: o=sp.kron(o,x,format='csr')
        cs.append(o)
    return cs
def L(n,d,z,t=1.0):
    c=ops(n); D=2**n; Id=sp.identity(D,format='csr')
    H=sum(-t*(c[j].T@c[(j+1)%n]+c[(j+1)%n].T@c[j]) for j in range(n))
    # row-major vec: vec(A r B)=kron(A,B.T)
    Lm=-1j*(sp.kron(H,Id)-sp.kron(Id,H.T))
    for j in range(n): Lm=Lm+sp.kron(c[j],c[j].conj())
    Lm=Lm-n*sp.identity(D*D)
    for xi in itertools.combinations(range(n),d):
        o=Id
        for j in xi: o=o@c[j].T
        Lm=Lm+z*sp.kron(o,o.conj())
    return Lm.tocsc()
if __name__=="__main__":
  n=int(sys.argv[1]); d=int(sys.argv[2])
  for z in (1e-4,1e-2):
      Lm=L(n,d,z)
      if Lm.shape[0]<=4096: ev=np.linalg.eigvals(Lm.toarray())
      else: ev=sla.eigs(Lm,k=30,sigma=0.0+0j,which='LM',return_eigenvectors=False)  # near -n shifted? use shift-invert at top
      ev=ev[np.argsort(-ev.real)]
      pb=np.sort(np.linalg.eigvals(Nmat(n,d,z)).real)[::-1]-n
      print(f'n={n} d={d} z={z}: full top Re:',np.round(ev.real[:6]+n,5),' | P-block top (shifted):',np.round(pb[:3]+n,5))
