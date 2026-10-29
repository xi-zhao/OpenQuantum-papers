import numpy as np, scipy.sparse.linalg as sla, json, sys
from ex import liou, obs_vec, ops
V=10.;Om=2.0
# Test B input: quasi-potential Phi(n) = -(1/N) ln P_ss(k) from the exact steady state (real Delta, independent of EP)
def Pk(N,D):
    L=liou(N,D,Om,V); tr,nv=obs_vec(N); Ne,Jx,B,idx,Dd=ops(N)
    w,v=sla.eigs(L,k=1,sigma=1e-7,tol=1e-14); r=v[:,0]/(tr@v[:,0])
    d=Ne.shape[0]; full=np.zeros(d*d,complex); full[idx]=r; R=full.reshape(d,d,order='F')
    ne=np.round(np.diag(Ne).real).astype(int); P=np.zeros(N+1)
    for i in range(d): P[ne[i]]+=R[i,i].real
    return P
N=int(sys.argv[1]); out={}
for D in [3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.62,3.65,3.7]:
    P=Pk(N,D); out[str(D)]=P.tolist(); print(N,D,'min P=%.2e'%P.min(),flush=True)
json.dump(out,open(f'../data/phi_N{N}.json','w'))
