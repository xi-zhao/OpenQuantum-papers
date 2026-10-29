# same base point & tangent velocity: big loop encircling EP vs small loop not encircling it
import numpy as np
from model import *
G=1.0;g=0.2;p0=2*np.pi/3
x0,y0=loop(0.0)            # base point (0.433,0.25)
cx,cy=0.0,0.5
for R,lab in [(0.5,'big (encircles EP)'),(0.1,'small (no EP)'),(0.02,'tiny')]:
    ccx=x0+(cx-x0)*R/0.5; ccy=y0+(cy-y0)*R/0.5
    a=R; b=R; c=ccy; # loop(phi)= a sin(phi+p0)+? need center x offset
    import model
    def lp(phi,R=R,ccx=ccx,ccy=ccy): return ccx+R*np.sin(phi+p0), ccy+R*np.cos(phi+p0)
    model.loop=lp
    enc=np.hypot(ccx-0,ccy-0.5)<R
    for T0 in [300,1000]:
        T=T0*R/0.5   # equal speed at base point
        c=tdist(evolve(T,1,G,g,nst=int(20*T)+2000),evolve(T,-1,G,g,nst=int(20*T)+2000))
        print(lab,'encircles=',enc,'T=%.1f'%T,'C=%.5e'%c,' C*T0=%.4f'%(c*T0))
