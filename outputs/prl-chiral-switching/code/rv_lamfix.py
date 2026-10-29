import numpy as np, json, scipy.sparse.linalg as sla
from ex import liou
V=10.; Om=2.0
out=[]
for N in [20,30,40,50,60,70]:
    for D in [3.0,3.2,3.3,3.4,3.5,3.55,3.6,3.65]:
        L=liou(N,D,Om,V); w=[]
        for sig in [-0.4,-0.8,-0.7+2.8j]:
            w+=list(sla.eigs(L,k=6,sigma=sig,return_eigenvectors=False,tol=1e-10))
        w=np.unique(np.round(np.array(w),7)); w=w[np.argsort(-w.real)]; w=w[1:6]
        out.append([N,D,[[x.real,x.imag] for x in w]]); print(N,D,np.round(w,4),flush=True)
        json.dump(out,open('../data/lamfix.json','w'))
