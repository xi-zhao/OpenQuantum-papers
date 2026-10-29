import numpy as np, scipy.sparse.linalg as sla, json
from ex import liou, ops
V=10.;Om=2.0
def diagpop(N,v):
    Ne,Jx,B,idx,D=ops(N); d=Ne.shape[0]; full=np.zeros(d*d,complex); full[idx]=v; R=full.reshape(d,d,order='F')
    # population distribution over excitation number k (Dicke basis grouped by Ne eigenvalue)
    ne=np.round(np.diag(Ne).real).astype(int); P=np.zeros(N+1)
    for i in range(d): P[ne[i]]+=R[i,i].real
    return P
out={}
for N in [40]:
    for D in [3.15,3.20,3.25,3.30,3.40]:
        L=liou(N,D,Om,V); w,v=sla.eigs(L,k=8,sigma=-0.6,tol=1e-12)
        sel=[i for i in np.argsort(-w.real) if abs(w[i].imag)<1e-7 and w[i].real<-0.05][:2]
        res=[]
        for i in sel:
            P=diagpop(N,v[:,i]); k=np.arange(N+1)/N; a=np.abs(P)/np.abs(P).sum()
            res.append(dict(lam=w[i].real,nbar=float((a*k).sum()),peak=float(k[np.argmax(abs(P))]),P=P.tolist()))
        out[f'{N},{D}']=res
        print(N,D,' | '.join('lam=%.4f  <n>_|P|=%.3f peak=%.3f'%(r['lam'],r['nbar'],r['peak']) for r in res),flush=True)
json.dump(out,open('../data/g1_modes.json','w'))
