import numpy as np, scipy.sparse.linalg as sla, sys, json
from ryd import PeriodicChain, lanczos_spectral
Dc,Vc=-1.5011,-1.6260
E7=np.array([1,1.28558,1.87939,1.96962,2.53209,2.87939,3.70167])
def run(L,th,t,m=400):
    pc=PeriodicChain(L); a=np.deg2rad(th); D=Dc+t*np.cos(a); V=Vc+t*np.sin(a)
    H0=pc.H(D,V,0); e,v=sla.eigsh(H0,k=1,which='SA',tol=1e-13); E0=e[0]; g=v[:,0]
    c0,_=pc.ops(0); o=c0*g; o-=(g@o)*g
    e0,w0=lanczos_spectral(H0,o,m)
    Hp=pc.H(D,V,1); cp,tgt=pc.ops(1); op=np.zeros(Hp.shape[0]); mk=tgt>=0
    np.add.at(op,tgt[mk],cp[mk]*g[mk]); ep,wp=lanczos_spectral(Hp,op,m)
    return e0-E0,w0/L,ep-E0,wp/L
L=int(sys.argv[1]); th=float(sys.argv[2]); t=float(sys.argv[3])
e0,w0,ep,wp=run(L,th,t)
m1=ep[np.argmax(wp>1e-8)]
res=dict(L=L,th=th,t=t,e0=e0.tolist(),w0=w0.tolist(),ep=ep.tolist(),wp=wp.tolist())
json.dump(res,open(f'../data/hispec_L{L}_th{th:.0f}_t{t}.json','w'))
for name,e,w in (('pi',ep,wp),('0',e0,w0)):
    s=(w>1e-7)&(e<4.2*m1)
    print(name, ' '.join(f"{x/m1:.3f}[{y:.1e}]" for x,y in zip(e[s],w[s])))
