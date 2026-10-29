import numpy as np, mpmath as mp
from pblock import Nmat
from math import comb
mp.mp.dps=50
# (a) zero root: algebraic r, geometric multiplicity = dim ker
for n,d in [(7,2),(9,3),(10,2),(11,3),(13,4)]:
    N=mp.matrix(Nmat(n,d,1.0).tolist())
    # rank via SVD in mp
    s=mp.svd_r(N,compute_uv=False)
    ker=sum(1 for x in s if x<mp.mpf(10)**-30)
    print(f'n={n} d={d} r={(n+1)%(d+1)} dimker={ker}')
