import numpy as np, itertools
from scipy.optimize import brentq, minimize
from theory2 import lnT, LN2
from exponent import E, qc
L=np.log2
# (1) extrapolate exponent: fit -log2T/n = a + b log2(n)/n + c/n
for r in [0.3,0.35,0.4,0.5,0.7]:
    ns=np.array([100,200,400,800]); y=np.array([-lnT(n,r,40001)/LN2/n for n in ns])
    A=np.stack([np.ones_like(ns,float),np.log2(ns)/ns,1/ns],1); a=np.linalg.lstsq(A,y,rcond=None)[0]
    print("extrap",r,round(a[0],4),"E",round(E(r),4),"b",round(a[1],3))
# (2) T at r_c
rc=4*qc/3
print("T(rc)",[round(np.exp(lnT(n,rc,40001)),4) for n in [25,50,100,200,400]])
# (3) general formula over 4-letter types for depolarizing
def Eg(p):
    def f(x):
        t=np.exp(x)/np.exp(x).sum(); H=-(t*L(t)).sum(); D=(t*L(t/p)).sum()
        return D+0.5*max(0,H-1)
    return min(minimize(f,np.random.randn(4),method='Nelder-Mead',options={'xatol':1e-10,'fatol':1e-12,'maxiter':20000}).fun for _ in range(8))
for r in [0.3,0.5,0.7]:
    q=3*r/4; print("simplex",r,round(Eg(np.array([1-q,q/3,q/3,q/3])),4),round(E(r),4))
