import mpmath as mp
from pblock import Nmat
mp.mp.dps=80
for d in (2,3):
    for n in (12,20,30,40):
        ev=mp.eig(mp.matrix(Nmat(n,d,1.0).tolist()))[0]
        re=sorted([mp.re(e) for e in ev],reverse=True)
        k1=re[0]; k2=next(x for x in re if x<k1-mp.mpf(10)**-20)
        evm=mp.eig(mp.matrix(Nmat(n-1,d,1.0).tolist()))[0]; km=max(mp.re(e) for e in evm)
        print(d,n,'kappa_n=%.5f 2nd Re in K_n=%.5f kappa_{n-1}=%.5f  internal/inter=%.4f (asym %.4f)'%(k1,k2,km,(k1-k2)/(k1-km),(2*(d+1))**0.5))
