# similarity-symmetrize: N is tridiagonal-like? not for d>=2. Use diagonal scaling S=diag(s_k) to balance, then double eig.
import numpy as np
from math import comb, lgamma
from pblock import Nmat
from scipy.linalg import eigvals
def mu_bal(n,d,z):
    N=Nmat(n,d,z)
    # balance with s_k chosen so super-diagonal entries become geometric mean scale: log s_k from WKB p(x)
    from scipy.linalg import matrix_balance
    B,_=matrix_balance(N,permute=False)
    return np.max(eigvals(B).real)
for d in (2,3):
    for zeta in (1e-4,1,100):
        print(d,zeta,[round(mu_bal(n,d,zeta/n**(d-1))/n,4) for n in (20,40,80,160,320)])
