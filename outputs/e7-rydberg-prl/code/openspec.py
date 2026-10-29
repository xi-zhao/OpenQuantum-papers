import numpy as np, scipy.sparse.linalg as sla, sys, json
from ryd import OpenChain, lanczos_spectral
Dc,Vc=-1.5011,-1.6260
L=int(sys.argv[1]); th=float(sys.argv[2]); t=float(sys.argv[3]); m=int(sys.argv[4]) if len(sys.argv)>4 else 300
oc=OpenChain(L); a=np.deg2rad(th); D=Dc+t*np.cos(a); V=Vc+t*np.sin(a)
H=oc.H(D,V); e,v=sla.eigsh(H,k=1,which='SA',tol=1e-13); E0=e[0]; g=v[:,0]
nb=oc.bits.astype(float)
dens=(nb*g[:,None]**2).sum(0)
out=dict(L=L,th=th,t=t,D=D,V=V,dens=dens.tolist())
for name,f in (('u',np.ones(L)),('s',(-1.0)**np.arange(L))):
    o=(nb@f)*g; o-=(g@o)*g
    ee,ww=lanczos_spectral(H,o,m)
    out[name]=dict(e=(ee-E0).tolist(),w=ww.tolist())
json.dump(out,open(f'../data/open_L{L}_th{th:.0f}_t{t}.json','w'))
for name in 'us':
    ee=np.array(out[name]['e']); ww=np.array(out[name]['w']); s=(ww>1e-4*ww.max())&(ee<4)
    print(name,' '.join(f"{x:.3f}[{y:.1e}]" for x,y in zip(ee[s],ww[s])))
print('dens',np.round(dens,3))
