import numpy as np
from arc import Rubidium87
A=Rubidium87(); n=53
for T in (0,300):
    print('T',T,'lifetime us', A.getStateLifetime(n,0,0.5,temperature=T,includeLevelsUpTo=n+40)*1e6)
tot=0
for nn,jj in [(52,1.5),(53,1.5),(52,0.5),(53,0.5),(51,1.5),(54,1.5)]:
    r0=A.getTransitionRate(n,0,0.5,nn,1,jj,temperature=0); r3=A.getTransitionRate(n,0,0.5,nn,1,jj,temperature=300)
    f=A.getTransitionFrequency(n,0,0.5,nn,1,jj)
    d=A.getReducedMatrixElementJ(n,0,0.5,nn,1,jj)  # in a0 e
    print(f'53S->{nn}P{jj}: f={f/1e9:.2f}GHz G0={r0:.2f} G300={r3:.0f}/s nth={r3/r0-1:.0f} <||d||>={d:.0f} ea0')
# C3 for |nS,n'P> exchange (radial, angle-averaged scale): C3 ~ d^2/(4 pi eps0) with d = reduced/sqrt(3) approx
a0=5.29177e-11; e=1.602e-19; eps0=8.854e-12; h=6.626e-34
for nn,jj in [(52,1.5),(53,1.5)]:
    d=A.getReducedMatrixElementJ(n,0,0.5,nn,1,jj)*a0*e/np.sqrt(3)
    C3=d**2/(4*np.pi*eps0)/h  # Hz m^3
    for r in (2e-6,5e-6,12e-6,50e-6,200e-6):
        print(f'  {nn}P{jj}: C3/r^3 at {r*1e6:.0f}um = {C3/r**3/1e6:.3f} MHz')
