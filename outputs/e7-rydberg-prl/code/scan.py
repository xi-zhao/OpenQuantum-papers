import numpy as np, scipy.sparse.linalg as sla, sys, json
from ryd import PeriodicChain, lanczos_spectral
Dc,Vc=-1.502,-1.627
L=int(sys.argv[1]); pc=PeriodicChain(L)
def spec(D,V,m=250):
    H0=pc.H(D,V,0); e,v=sla.eigsh(H0,k=1,which='SA',tol=1e-13); E0=e[0]; g=v[:,0]
    c0,_=pc.ops(0); o=c0*g; o-= (g@o)*g
    e0,w0=lanczos_spectral(H0,o,m)
    Hp=pc.H(D,V,1); cp,tgt=pc.ops(1); op=np.zeros(Hp.shape[0]); mk=tgt>=0
    np.add.at(op,tgt[mk],cp[mk]*g[mk])
    ep,wp=lanczos_spectral(Hp,op,m)
    return e0-E0,w0,ep-E0,wp
for th in [float(x) for x in sys.argv[2].split(',')]:
  for t in [0.15,0.3]:
    ang=np.deg2rad(th); D=Dc+t*np.cos(ang); V=Vc+t*np.sin(ang)
    e0,w0,ep,wp=spec(D,V)
    sel0=w0>1e-3*w0.max(); selp=wp>1e-4*wp.max()
    m1=ep[selp][0]
    print(f"th={th:5.0f} t={t:.2f} D={D:.3f} V={V:.3f} m1={m1:.4f}")
    print("  pi:",' '.join(f"{x/m1:.3f}({w/wp.max():.1e})" for x,w in zip(ep[selp][:7],wp[selp][:7])))
    print("  0 :",' '.join(f"{x/m1:.3f}({w/w0.max():.1e})" for x,w in zip(e0[sel0][:7],w0[sel0][:7])),flush=True)
