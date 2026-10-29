import numpy as np
from mf import roots
V=10.
def ell(th,c=(4.9,3.3),a=(1.5,1.7),tilt=0.0):
    return c[0]+a[0]*np.cos(th), c[1]+a[1]*np.sin(th)
th=np.linspace(0,2*np.pi,721)
nb=[len(roots(*ell(t),V)) for t in th]
s=''.join('3' if x==3 else '.' for x in nb); print(s[::6])
