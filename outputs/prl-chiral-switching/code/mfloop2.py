import numpy as np
from scipy.integrate import solve_ivp
from mfdyn import rhs, steady
from mf import roots
V=10.
p0=(4.1,2.0)
def run(c0,R,T,sgn,br):
    R=np.array(R); th0=np.arctan2((p0[1]-c0[1])/R[1],(p0[0]-c0[0])/R[0])
    assert abs(c0[0]+R[0]*np.cos(th0)-p0[0])<1e-9 and abs(c0[1]+R[1]*np.sin(th0)-p0[1])<1e-9
    path=lambda t:(c0[0]+R[0]*np.cos(th0+sgn*2*np.pi*t/T), c0[1]+R[1]*np.sin(th0+sgn*2*np.pi*t/T), V)
    s=solve_ivp(rhs,[0,T],steady(*p0,V,br),args=(path,),rtol=1e-8,atol=1e-10,method='LSODA')
    return s.y[0,-1]
# loop A: encircles cusp. loop B: horizontal hysteresis loop at low Om, not enclosing cusp
for name,c0,R in [('encircle',(4.9,3.4),(1.6,1.6)),('nocusp',(4.1,1.6),(1.3,0.4*np.sqrt(1)) )]:
    R=np.array(R)
    # rescale R[1] so loop passes through p0
    if name=='nocusp': R=np.array([1.3, 0.4]); 
    for br in [0,2]:
        for sgn in [1,-1]:
            try: print(name,br,sgn,round(run(c0,R,500,sgn,br),4))
            except AssertionError: print('loop misses p0'); break
