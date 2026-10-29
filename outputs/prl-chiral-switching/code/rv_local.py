import numpy as np, json
from mf import roots
from mfdyn import rhs, steady
V=10.; Om=2.0
# MF slow-manifold (centre-manifold) eigenvalue on the high branch approaching the left fold, vs exact lambda's
Df=3.59963
def Jev(D,br):
    y=np.array(steady(D,Om,V,br)); f=lambda y: np.array(rhs(0,y,lambda t:(D,Om,V))); h=1e-7
    J=np.array([(f(y+h*e)-f(y-h*e))/(2*h) for e in np.eye(3)]).T; return np.sort_complex(np.linalg.eigvals(J))
for D in [3.0,3.2,3.4,3.55,3.59,3.6,3.62,3.65,3.7]:
    r=roots(D,Om,V); s='  '.join(f"br{k}(n={r[k]:.3f}):{np.round(Jev(D,k),3)}" for k in range(len(r)))
    print(D,s)
