import numpy as np
from lvl import roots3
from matplotlib.path import Path
xs=np.concatenate([np.linspace(0,0.02,200),np.linspace(0.02,0.3,300)[1:]])
# cusp: smallest Om with bistability
def has(O):
    return any(len(roots3(D,O,xs=xs))>=3 for D in np.linspace(-33,-27,121))
lo,hi=6.5,7.6
for _ in range(10):
    m=(lo+hi)/2; (hi:=m) if has(m) else (lo:=m)
print('cusp Om_p ~',round(hi,3))
Ds=np.linspace(-33,-27,601); b=[D for D in Ds if len(roots3(D,hi+0.02,xs=xs))>=3]; print('cusp Delta ~',np.mean(b))
Om0,pp,D0,dl=10.58,0.31,-52.80,14.92
th=np.linspace(0,2*np.pi,400)
loop=Path(np.c_[D0+dl*np.sin(th),Om0*np.sqrt(1+pp*np.cos(th))])
print('cusp inside experimental loop?',loop.contains_point((np.mean(b),hi)))
print('loop Om range',Om0*np.sqrt(1-pp),Om0*np.sqrt(1+pp),'Delta range',D0-dl,D0+dl)
