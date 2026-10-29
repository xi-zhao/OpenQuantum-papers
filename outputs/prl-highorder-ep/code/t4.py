import mpmath as mp
from pblock import Nmat
mp.mp.dps=120
for d,n,zeta in [(2,80,100),(2,160,100),(2,40,100),(3,160,1e-4),(2,160,1e-4)]:
    z=zeta/n**(d-1)
    ev=mp.eig(mp.matrix(Nmat(n,d,z).tolist()))[0]
    print(d,n,zeta,'mu/n=%.4f'%(float(max(mp.re(e) for e in ev))/n))
