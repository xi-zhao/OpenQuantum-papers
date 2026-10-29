import sympy as sp, numpy as np
from math import comb
from pblock import Nmat
from scipy.linalg import eigvals
mu_=sp.Symbol('mu'); z=sp.Symbol('z',positive=True)
# (a) exact char poly structure and nonzero lowest coefficient
bad=[]
for d in range(1,6):
    for n in range(d,25):
        M=sp.zeros(n+1,n+1)
        for k in range(1,n+1): M[k-1,k]=n-k+1
        for k in range(0,n+1-d): M[k+d,k]=z*comb(k+d,d)
        p=sp.Poly((mu_*sp.eye(n+1)-M).det(method='berkowitz').expand(),mu_)
        r=(n+1)%(d+1)
        lowdeg=min(m[0] for m in p.monoms())
        if lowdeg!=r: bad.append((n,d,lowdeg,r))
print('mismatches (n,d,lowdeg,r):',bad)
