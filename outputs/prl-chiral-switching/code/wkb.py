import numpy as np, json
from scipy.integrate import quad
from scipy.optimize import brentq
V=10.;Om=2.0
# WKB action of the rate-equation (birth-death) reduction with the SAME fixed points as the Bloch MF:
#   w+(n)=(1-n) r(n), w-(n)=n (1+r(n)), r(n)=(Om^2/4)/((D-Vn)^2+1/4)
# quasi-potential Phi'(n)=ln(w-/w+);  barrier between stable branch and saddle
def r(n,D): return (Om**2/4)/((D-V*n)**2+0.25)
def dphi(n,D): return np.log(n*(1+r(n,D))/((1-n)*r(n,D)))
def fps(D):
    g=lambda n: (1-n)*r(n,D)-n*(1+r(n,D)); xs=np.linspace(1e-4,0.499,3000); f=g(xs)
    return [brentq(g,xs[i],xs[i+1]) for i in np.flatnonzero(np.sign(f[:-1])!=np.sign(f[1:]))]
for D in [3.62,3.7,3.8,3.9,4.0]:
    f=fps(D)
    if len(f)==3:
        lo,sad,hi=f
        print(D,'hi->lo barrier %.4f   lo->hi barrier %.4f'%(quad(dphi,hi,sad,args=(D,))[0]*-1,quad(dphi,lo,sad,args=(D,))[0]))
