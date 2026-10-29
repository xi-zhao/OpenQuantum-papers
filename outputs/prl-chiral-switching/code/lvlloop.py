import numpy as np
from lvl import roots3
# experimental loop of Xie et al. (Omega_p, Delta_c) plane, Methods: Om0=10.58, p=0.31, D0=-52.80, delta=14.92 (2pi MHz)
Om0,p,D0,dl=10.58,0.31,-52.80,14.92
th=np.linspace(0,2*np.pi,241)
pts=[(D0+dl*np.sin(t),Om0*np.sqrt(1+p*np.cos(t))) for t in th]
nb=[len(roots3(D,O)) for D,O in pts]
print(''.join('3' if k>=3 else '.' for k in nb))
e=np.flatnonzero(np.diff(np.array(nb)>=3)); print('crossings theta=',np.round(th[e],3))
for i in e: print(np.round(pts[i],2))
