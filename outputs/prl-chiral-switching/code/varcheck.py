import numpy as np, scipy.sparse.linalg as sla
from ex import liou, obs_vec, ops
V=10.
def exact_var(N,D,Om):
    L=liou(N,D,Om,V); tr,nv=obs_vec(N); Ne,Jx,B,idx,Dd=ops(N)
    w,v=sla.eigs(L,k=1,sigma=1e-7); r=v[:,0]/(tr@v[:,0])
    n2=(Ne@Ne).T.flatten(order='F')[idx]; m=(nv@r).real; return m/N, ((n2@r).real-m*m)/N
for D,Om in [(4.9,2.0),(3.0,2.0),(2.0,1.0)]:
    print(D,Om,[ (N,)+tuple(np.round(exact_var(N,D,Om),4)) for N in (20,40,60)])
