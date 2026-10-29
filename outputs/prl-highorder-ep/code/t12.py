# Is the second-largest Re of the whole spectrum mu_{n-1}? i.e. non-symmetric S_n sectors of G_n lie below mu_{n-1}
import numpy as np, itertools
from t11 import G
from pblock import Nmat
from scipy.linalg import eigvals
mu=lambda m,d: np.max(eigvals(Nmat(m,d,1.0)).real) if m>=0 else -np.inf
for d in (1,2,3,4):
    for n in range(d+1,11):
        ev=np.sort(eigvals(G(n,d,1.0)).real)[::-1]
        sym=np.sort(eigvals(Nmat(n,d,1.0)).real)[::-1]
        # remove symmetric-sector eigenvalues once each
        rest=list(ev)
        for s in sym:
            i=int(np.argmin(np.abs(np.array(rest)-s))); rest.pop(i)
        print(f'd={d} n={n}: mu_n={sym[0]:.4f} mu_n-1={mu(n-1,d):.4f} max nonsym(G_n)={max(rest):.4f}  2nd sym={sym[1]:.4f}')
