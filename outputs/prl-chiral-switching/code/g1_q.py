import numpy as np, scipy.sparse.linalg as sla, json
from ex import liou, ops
V=10.;Om=2.0
# compare the "ghost" mode's population profile with the MF slow flow: does its weight sit in the ghost region of the
# vanished low branch (n ~ n_fold_low ~ 0.155 at the fold)?
def pop(N,v):
    Ne,Jx,B,idx,D=ops(N); d=Ne.shape[0]; full=np.zeros(d*d,complex); full[idx]=v; R=full.reshape(d,d,order='F')
    ne=np.round(np.diag(Ne).real).astype(int); P=np.zeros(N+1)
    for i in range(d): P[ne[i]]+=R[i,i].real
    return P
for N in [40,60]:
    for D in [3.3,3.4,3.5]:
        L=liou(N,D,Om,V); w,v=sla.eigs(L,k=8,sigma=-0.3,tol=1e-12)
        i=[j for j in np.argsort(-w.real) if abs(w[j].imag)<1e-7 and w[j].real<-1e-4][0]
        P=pop(N,v[:,i]); P/=np.max(abs(P)); k=np.arange(N+1)/N
        mx=k[np.argmax(P)]; mn=k[np.argmin(P)]
        print(N,D,'lam=%.4f  max at n=%.3f  min at n=%.3f'%(w[i].real,mx,mn), ' P(n) coarse:',np.round(P[::max(1,N//10)],2))
