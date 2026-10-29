import numpy as np, json, scipy.sparse.linalg as sla
from multiprocessing import Pool
from ex import liou, obs_vec, ops
V=10.;Om=2.0
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
nf=json.load(open('../data/normalform.json'))['nf']
def Pk(N,D):
    L=liou(N,D,Om,V); tr,nv=obs_vec(N); Ne,Jx,B,idx,Dd=ops(N)
    w,v=sla.eigs(L,k=1,sigma=1e-9,tol=1e-15); r=v[:,0]/(tr@v[:,0])
    d=Ne.shape[0]; full=np.zeros(d*d,complex); full[idx]=r; R=full.reshape(d,d,order='F')
    ne=np.round(np.diag(Ne).real).astype(int); P=np.zeros(N+1)
    for i in range(d): P[ne[i]]+=R[i,i].real
    return P
def job(N):
    D=EP[N][1]; P=Pk(N,D); return [N,D,P.tolist()]
if __name__=="__main__":
    with Pool(6) as Pp: res=Pp.map(job,[20,24,28,32,36,40,44,48])
    json.dump(res,open('../data/phiEP.json','w')); print('done')
