import numpy as np, json, scipy.sparse.linalg as sla, sys
from ex import liou
V=10.; Om=2.0; N=int(sys.argv[1]); out=[]
for D in np.linspace(2.9,3.7,33):
    L=liou(N,D,Om,V); w=[]
    for sig in [-0.3,-0.75,-0.75+2.8j,-1.1]:
        w+=list(sla.eigs(L,k=6,sigma=sig,return_eigenvectors=False,tol=1e-10))
    w=np.unique(np.round(np.array(w),7)); w=w[w.real<-1e-8]; w=w[np.argsort(-w.real)][:6]
    out.append([D,[[x.real,x.imag] for x in w]]); print(N,round(D,3),np.round(w[:4],4),flush=True)
json.dump(out,open(f'../data/dense_N{N}.json','w'))
