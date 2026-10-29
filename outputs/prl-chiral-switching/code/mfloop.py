import numpy as np
from scipy.integrate import solve_ivp
from mfdyn import rhs, steady
from mf import roots
V=10.
c0=(4.9,3.4); R=1.6
p0=(4.1,2.0); th0=np.arctan2(p0[1]-c0[1],p0[0]-c0[0])
def mkpath(T,sgn):
    def path(t):
        th=th0+sgn*2*np.pi*t/T
        return c0[0]+R*np.cos(th), c0[1]+R*np.sin(th), V
    return path
print('start roots',roots(*p0,V))
for T in [50,200,1000]:
  for br,name in [(0,'low'),(2,'high')]:
    for sgn,dn in [(1,'CCW'),(-1,'CW')]:
        y0=steady(p0[0],p0[1],V,br); path=mkpath(T,sgn)
        s=solve_ivp(rhs,[0,T],y0,args=(path,),rtol=1e-8,atol=1e-10,method='LSODA')
        print(T,name,dn,round(y0[0],4),'->',round(s.y[0,-1],4))
