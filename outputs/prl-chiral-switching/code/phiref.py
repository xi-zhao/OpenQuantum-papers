import numpy as np, json, scipy.sparse.linalg as sla
from ex import liou, obs_vec, ops
V=10.;Om=2.0
def Pk(N,D):
    L=liou(N,D,Om,V); tr,nv=obs_vec(N); Ne,Jx,B,idx,Dd=ops(N)
    w,v=sla.eigs(L,k=1,sigma=1e-9,tol=1e-15); r=v[:,0]/(tr@v[:,0])
    d=Ne.shape[0]; full=np.zeros(d*d,complex); full[idx]=r; R=full.reshape(d,d,order='F')
    ne=np.round(np.diag(Ne).real).astype(int); P=np.zeros(N+1)
    for j in range(d): P[ne[j]]+=R[j,j].real
    return P
out={}
for D,ns in [(3.25,[0.15693]),(3.7,[0.25,0.30,0.35])]:
    vals={N:Pk(N,D) for N in (30,40,50)}
    for n in ns:
        ph=[]
        for N in (30,40,50):
            p=vals[N]; lp=np.log(np.clip(p,1e-300,None)); k=np.arange(N+1)/N; ph.append(-(np.interp(n,k,lp)-lp.max())/N)
        ext=np.polyfit([1/30,1/40,1/50],ph,1)[1]; out[f'{D},{n}']=dict(PhiN=ph,ext=ext); print(D,n,np.round(ph,4),'ext=%.4f'%ext,flush=True)
json.dump(out,open('../data/phiref_gateA.json','w'))
