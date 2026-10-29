import numpy as np, scipy.sparse.linalg as sla, json
from ryd import PeriodicChain
from scipy.optimize import fsolve
Ls=(16,20,24,28,32)
ch={L:PeriodicChain(L) for L in Ls}
def lev(L,D,V,k0=4,kp=3):
    pc=ch[L]
    e0=np.sort(sla.eigsh(pc.H(D,V,0),k=k0,which='SA',tol=1e-12,return_eigenvectors=False))
    ep=np.sort(sla.eigsh(pc.H(D,V,1),k=kp,which='SA',tol=1e-12,return_eigenvectors=False))
    return e0-e0[0], ep-e0[0]
def F(p,L1,L2):
    D,V=p
    a0,ap=lev(L1,D,V,2,1); b0,bp=lev(L2,D,V,2,1)
    return [L1*ap[0]-L2*bp[0], L1*a0[1]-L2*b0[1]]
res=[]
guess=(-1.5,-1.63)
for L1,L2 in zip(Ls[:-1],Ls[1:]):
    sol=fsolve(F,guess,args=(L1,L2),xtol=1e-10); guess=sol
    e0,ep=lev(L2,*sol,6,4)
    xs=ep[0]*L2  # 2 pi v x_sigma
    v=xs/(2*np.pi*3/40)
    print(f"{L1}-{L2}: Dc={sol[0]:.6f} V2c={sol[1]:.6f} v={v:.4f}  k0 x:",np.round(e0[1:]*L2/(2*np.pi*v),4)," kpi x:",np.round(ep*L2/(2*np.pi*v),4),flush=True)
    res.append(dict(L1=L1,L2=L2,Dc=sol[0],V2c=sol[1],v=v))
json.dump(res,open('../data/tci.json','w'),indent=1)
