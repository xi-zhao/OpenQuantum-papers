# effective number of independent Gaussian channels for BBR transitions from Rb 53S1/2 m=+1/2
import numpy as np, json
from arc import Rubidium87
from sympy.physics.quantum.cg import CG
from sympy import S
A=Rubidium87(); n=53
T=300
G0=1/A.getStateLifetime(n,0,0.5,temperature=0,includeLevelsUpTo=n+40)
G3=1/A.getStateLifetime(n,0,0.5,temperature=T,includeLevelsUpTo=n+40)
ch=[]
for nn in range(5,n+40):
  for jj in (0.5,1.5):
    try:
      r0=A.getTransitionRate(n,0,0.5,nn,1,jj,temperature=0); r3=A.getTransitionRate(n,0,0.5,nn,1,jj,temperature=T)
    except Exception: continue
    f=A.getTransitionFrequency(n,0,0.5,nn,1,jj)
    bbr=r3-r0
    if bbr<=0: continue
    # split by polarization q using CG weights for m=1/2 (rates for isotropic field: total over q weighted by |CG|^2 normalized)
    w=[float(CG(S(1)/2,S(1)/2,1,q,S(jj),S(1)/2+q).doit())**2 if abs(S(1)/2+q)<=jj else 0 for q in (-1,0,1)]
    w=np.array(w)/sum(w)
    for q,wq in zip((-1,0,1),w):
        if wq>0: ch.append(dict(nn=nn,j=jj,q=q,f=f,rate=bbr*wq))
r=np.array([c['rate'] for c in ch])
nu=r.sum()**2/(r**2).sum()
print('Gamma_sp(0K)=%.0f /s  Gamma_BBR=%.0f /s  sum channels=%.0f'%(G0,G3-G0,r.sum()))
print('nu_eff=%.2f'%nu)
top=sorted(ch,key=lambda c:-c['rate'])[:8]
for c in top: print(c)
json.dump(dict(G_sp=G0,G_bbr=G3-G0,nu_eff=nu,channels=ch),open('../data_nu.json','w'),default=float)
