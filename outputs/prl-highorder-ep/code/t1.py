import numpy as np, mpmath as mp
from pblock import *
mp.mp.dps=60
z=1e-5
# (1) zero roots & homogeneity, exact in high precision
for n,d in [(4,1),(4,2),(4,3),(4,4),(7,2),(9,3),(10,2)]:
    N=mp.matrix(Nmat(n,d,1.0).tolist())
    Z=mp.matrix(Nmat(n,d,z).tolist())
    ev=mp.eig(Z)[0]
    small=sum(1 for e in ev if abs(e)<mp.mpf(10)**-25)
    em=max(ev,key=lambda e: mp.re(e))
    ev1=mp.eig(N)[0]; e1=max(ev1,key=lambda e: mp.re(e))
    print(f'n={n} d={d}: #zero roots={small} (pred r={(n+1)%(d+1)})  mu_max(z)/z^(1/(d+1))={float(mp.re(em))/z**(1/(d+1)):.6f}  mu_max(z=1)={float(mp.re(e1)):.6f}')
