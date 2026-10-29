import numpy as np, scipy.sparse.linalg as sla, json, sys
from ex import liou, obs_vec, ops
V=10.;Om=2.0
# Ground truth for the WKB action: the exact left eigenvector of the switching mode / exact P_ss(k)
# The steady-state population distribution P_ss(k) = exp(-N Phi(k/N)) defines the quasi-potential, whose
# value at the saddle relative to the hi minimum IS the WKB instanton action (Freidlin–Wentzell / Doi-Peliti).
def Pk(N,D):
    L=liou(N,D,Om,V); tr,nv=obs_vec(N); Ne,Jx,B,idx,Dd=ops(N)
    w,v=sla.eigs(L,k=1,sigma=1e-9,tol=1e-15); r=v[:,0]/(tr@v[:,0])
    d=Ne.shape[0]; full=np.zeros(d*d,complex); full[idx]=r; R=full.reshape(d,d,order='F')
    ne=np.round(np.diag(Ne).real).astype(int); P=np.zeros(N+1)
    for i in range(d): P[ne[i]]+=R[i,i].real
    return P
N=int(sys.argv[1]); out={}
for D in [3.7,3.8,3.9,4.0,4.1]:
    P=Pk(N,D); out[str(D)]=P.tolist(); print(N,D,'minP=%.2e'%P.min(),flush=True)
json.dump(out,open(f'../data/Pk_N{N}.json','w'))
