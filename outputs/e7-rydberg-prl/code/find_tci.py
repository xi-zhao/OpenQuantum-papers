import numpy as np, scipy.sparse.linalg as sla
from ryd import PeriodicChain
from scipy.optimize import brentq
chains={L:PeriodicChain(L) for L in (20,24,28)}
def gaps(L,D,V):
    pc=chains[L]
    e0=np.sort(sla.eigsh(pc.H(D,V,0),k=2,which='SA',tol=1e-11,return_eigenvectors=False))
    ep=np.sort(sla.eigsh(pc.H(D,V,1),k=1,which='SA',tol=1e-11,return_eigenvectors=False))
    return ep[0]-e0[0], e0[1]-e0[0]
for V in [-1.2,-1.3,-1.4,-1.5,-1.6,-1.7,-1.8]:
    f=lambda D: 20*gaps(20,D,V)[0]-24*gaps(24,D,V)[0]
    Ds=np.linspace(-1.5,2.0,15); fv=[f(D) for D in Ds]
    s=[i for i in range(len(Ds)-1) if np.sign(fv[i])!=np.sign(fv[i+1])]
    for i in s:
        try: Dc=brentq(f,Ds[i],Ds[i+1],xtol=1e-7)
        except ValueError: continue
        out=[(L,)+gaps(L,Dc,V) for L in (20,24,28)]
        print(f"V2={V:5.2f} Dc={Dc:.5f}", " ".join(f"L{L}: R={g0/gp:.3f} Lgp={L*gp:.3f}" for L,gp,g0 in out),flush=True)
