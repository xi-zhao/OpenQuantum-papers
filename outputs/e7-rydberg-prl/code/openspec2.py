import numpy as np, scipy.sparse.linalg as sla, sys, json
from ryd import OpenChain, lanczos_spectral
Dc,Vc=-1.5011,-1.6260
L=int(sys.argv[1]); th=float(sys.argv[2]); t=float(sys.argv[3]); m=int(sys.argv[4])
oc=OpenChain(L); a=np.deg2rad(th); D=Dc+t*np.cos(a); V=Vc+t*np.sin(a)
H=oc.H(D,V); e,v=sla.eigsh(H,k=1,which='SA',tol=1e-13); E0=e[0]; g=v[:,0]
nb=oc.bits.astype(float); j=np.arange(L)
env={'flat':np.ones(L),'taper':np.sin(np.pi*(j+1)/(L+1))**2}
out=dict(L=L,th=th,t=t,D=D,V=V)
for en,w in env.items():
  for name,f in (('u',w),('s',w*(-1.0)**j)):
    o=(nb@f)*g; o-=(g@o)*g
    ee,ww=lanczos_spectral(H,o,m)
    out[f'{name}_{en}']=dict(e=(ee-E0).tolist(),w=(ww/np.sum(f**2)).tolist())
    s=(ww>1e-3*ww.max())&(ee-E0<3)
    print(f'{name}_{en}',' '.join(f"{x:.3f}[{y:.1e}]" for x,y in zip((ee-E0)[s],(ww/np.sum(f**2))[s])),flush=True)
json.dump(out,open(f'../data/open2_L{L}_th{th:.0f}_t{t}.json','w'))
