import numpy as np, scipy.sparse.linalg as sla, sys, json
from ryd import PeriodicChain
Dc,Vc=-1.5015,-1.6262
E7=dict(m2=2*np.cos(5*np.pi/18),m3=2*np.cos(np.pi/9),m4=2*np.cos(np.pi/18))
def levels(pc,D,V,n0=6,npi=4):
    H0=pc.H(D,V,0); e0,v0=sla.eigsh(H0,k=n0,which='SA',tol=1e-12); o=np.argsort(e0); e0,v0=e0[o],v0[:,o]
    g=v0[:,0]; c0,_=pc.ops(0)
    w0=(v0[:,1:].T@(c0*g))**2
    Hp=pc.H(D,V,1); ep,vp=sla.eigsh(Hp,k=npi,which='SA',tol=1e-12); o=np.argsort(ep); ep,vp=ep[o],vp[:,o]
    cp,tgt=pc.ops(1); op=np.zeros(Hp.shape[0]); mk=tgt>=0; np.add.at(op,tgt[mk],cp[mk]*g[mk])
    wp=(vp.T@op)**2
    return e0[1:]-e0[0],w0/pc.L,ep-e0[0],wp/pc.L
if __name__=="__main__":
    Ls=[int(x) for x in sys.argv[1].split(',')]
    pts=[(float(a),float(b)) for a,b in (p.split(':') for p in sys.argv[2].split(','))]
    out=[]
    for L in Ls:
        pc=PeriodicChain(L)
        for th,t in pts:
            ang=np.deg2rad(th); D=Dc+t*np.cos(ang); V=Vc+t*np.sin(ang)
            e0,w0,ep,wp=levels(pc,D,V)
            m1=ep[0]
            print(f"L={L} th={th:.0f} t={t:.3f} m1={m1:.5f} | pi:",' '.join(f"{x/m1:.4f}[{w:.2e}]" for x,w in zip(ep,wp)),
                  "| 0:",' '.join(f"{x/m1:.4f}[{w:.2e}]" for x,w in zip(e0,w0)),flush=True)
            out.append(dict(L=L,th=th,t=t,e0=e0.tolist(),w0=w0.tolist(),ep=ep.tolist(),wp=wp.tolist()))
    json.dump(out,open(sys.argv[3],'w'))
