import numpy as np
from scipy.integrate import solve_ivp
from mfdyn import rhs, steady
from mf import roots, cusp
from matplotlib.path import Path
V=10.; P0=(3.7038915919767983,1.6194358071727704)
Dc,Oc=cusp(V)[0][1],cusp(V)[0][2]
# loop D: start P0, go up inside the wedge, around the cusp at large radius, but come back and re-enter through the SAME (left) fold:
# parametrize as polygon: P0 -> exit left fold to the left (low Delta, high branch only) -> up and over the cusp -> down on the right side
# staying OUTSIDE the wedge -> cross back to left of the wedge below P0 ... must cross right fold if passing below wedge? test numerically
pts=[P0,(2.0,1.9),(2.0,5.5),(7.0,5.5),(7.0,2.5),(6.0,0.5),(2.6,0.6),(2.6,1.45),P0]
P=np.array(pts)
seg=np.cumsum(np.r_[0,np.hypot(*np.diff(P,axis=0).T)]); seg/=seg[-1]
def mk(T,s):
    def path(t):
        u=(t/T) if s>0 else 1-t/T
        return float(np.interp(u,seg,P[:,0])),float(np.interp(u,seg,P[:,1])),V
    return path
print('encloses cusp:',Path(P).contains_point((Dc,Oc)))
u=np.linspace(0,1,4001); nb=np.array([len(roots(*mk(1,1)(x)[:2],V)) for x in u])
e=np.flatnonzero(np.diff(nb)); print('fold crossings at u=',np.round(u[e],3),[np.round(mk(1,1)(u[i])[:2],2) for i in e])
for s in (1,-1):
    for br in (0,2):
        S=solve_ivp(rhs,[0,2000],steady(*P0,V,br),args=(mk(2000,s),),rtol=1e-8,atol=1e-10,method='LSODA')
        print('dir',s,'init',br,'->',round(S.y[0,-1],4))
