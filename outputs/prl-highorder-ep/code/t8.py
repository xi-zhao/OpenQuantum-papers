import numpy as np, scipy.sparse.linalg as sla, sys
from full import L as Lfull
from pblock import Nmat
from scipy.linalg import eigvals
def mu(n,d,z): return np.max(eigvals(Nmat(n,d,z)).real)
for n,d in [(4,2),(5,2),(6,2),(7,2),(5,3),(6,3),(7,3)]:
    z=1e-3
    Lm=Lfull(n,d,z)
    D=Lm.shape[0]
    if D<=4096:
        ev=np.linalg.eigvals(Lm.toarray())
    else:
        ev=sla.eigs(Lm,k=40,sigma=-n+mu(n,d,z)+1e-3,which='LM',return_eigenvectors=False,tol=1e-12)
    ev=ev[np.argsort(-ev.real)]+n
    top=ev.real[0]; sec=ev.real[ev.real<top-1e-8][0]
    print(f'n={n} d={d}: top={top:.6f} (mu_n={mu(n,d,z):.6f})  2nd={sec:.6f} (mu_n-1={mu(n-1,d,z):.6f})  gap={top-sec:.6f}')
