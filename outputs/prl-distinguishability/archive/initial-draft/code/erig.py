import numpy as np
from scipy.optimize import brentq, minimize_scalar
from exponent import E, qc
L=np.log2
h=lambda f: -f*L(f)-(1-f)*L(1-f)
D=lambda a,b: a*L(a/b)+(1-a)*L((1-a)/(1-b))
def Erig(r):
    q=3*r/4
    s=lambda f: -(1-f)*L(1-q)-f*L(q/3)
    g=lambda f: 2*s(f)-h(f)-f*L(3)        # exponent of sum p^2 over type f
    fs=(q*q/3)/((1-q)**2+q*q/3)
    def Gexp(ft):  # min_{f>=ft} g
        return g(max(ft,fs))
    def val(ft):
        b = D(ft,q) if ft<q else 0.0      # B = {f<ft}
        return min(b,(Gexp(ft)-1)/2)
    res=minimize_scalar(lambda ft:-val(ft),bounds=(1e-9,q),method='bounded',options={'xatol':1e-12})
    return -res.fun
for r in [0.26,0.3,0.4,0.5,0.6,0.62,0.7,0.9]:
    print(r,round(Erig(r),5),round(E(r),5))
