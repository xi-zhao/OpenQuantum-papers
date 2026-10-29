import numpy as np, scipy.sparse.linalg as sla, sys, json
from ryd import PeriodicChain, lanczos_spectral
Dc,Vc=-1.5011,-1.6260
th,t=float(sys.argv[1]),float(sys.argv[2]); a=np.deg2rad(th); D=Dc+t*np.cos(a); V=Vc+t*np.sin(a)
out=[]
for L in range(18,37,2):
    pc=PeriodicChain(L); m=min(400,max(150,pc.B.reps.size//4))
    H0=pc.H(D,V,0); e,v=sla.eigsh(H0,k=1,which='SA',tol=1e-13); E0=e[0]; g=v[:,0]
    c0,_=pc.ops(0); o=c0*g; o-=(g@o)*g; e0,w0=lanczos_spectral(H0,o,m)
    Hp=pc.H(D,V,1); cp,tgt=pc.ops(1); op=np.zeros(Hp.shape[0]); mk=tgt>=0
    np.add.at(op,tgt[mk],cp[mk]*g[mk]); ep,wp=lanczos_spectral(Hp,op,m)
    out.append(dict(L=L,e0=(e0-E0).tolist(),w0=(w0/L).tolist(),ep=(ep-E0).tolist(),wp=(wp/L).tolist()))
    print(L,flush=True)
json.dump(out,open(f'../data/lscan_th{th:.0f}_t{t}.json','w'))
