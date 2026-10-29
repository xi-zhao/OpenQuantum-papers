import numpy as np, mpmath as mp
from pblock import *
mp.mp.dps=80
def lead(n,d):
    ev=mp.eig(mp.matrix(Nmat(n,d,1.0).tolist()))[0]
    ev=sorted(ev,key=lambda e:-mp.re(e)); return float(mp.re(ev[0])), float(mp.re(ev[0])-mp.re(ev[1]))
for d in (1,2,3):
    print('d=',d)
    for n in (8,16,32,64,96):
        m,g=lead(n,d); Hm,gs=semiclassical(n,d,1.0)
        print('  n=%3d mu_max=%.4f  WKB=%.4f ratio=%.4f | gap=%.4f WKB=%.4f'%(n,m,Hm,m/Hm,g,gs))
