# conjecture: spec G_n(z) = U_j [C(n,j)-C(n,j-1)] x spec N_{n-2j}(z)
import numpy as np
from math import comb
from t11 import G
from pblock import Nmat
from scipy.linalg import eigvals
def key(v): return np.sort_complex(np.round(v,7))
for d in (1,2,3):
    for n in range(d,10):
        ev=eigvals(G(n,d,0.37))
        pred=[]
        for j in range(0,n//2+1):
            mult=comb(n,j)-(comb(n,j-1) if j>0 else 0)
            m=n-2*j
            pred+=list(eigvals(Nmat(m,d,0.37)) if m>=d else eigvals(Nmat(m,d,0.0)))*mult
        a,b=key(ev),key(np.array(pred))
        ok=len(a)==len(b) and np.max(np.abs(a-b))<1e-5
        print(d,n,len(a),len(b),'OK' if ok else 'FAIL %.2e'%np.max(np.abs(a-b)) if len(a)==len(b) else 'len')
