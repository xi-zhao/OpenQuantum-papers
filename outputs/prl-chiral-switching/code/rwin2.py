import numpy as np, scipy.sparse.linalg as sla, json, sys
from multiprocessing import Pool
from ex import liou
V=10.;Om=2.0
c0={28:4.9966,30:4.987,32:4.9785,34:4.970,36:4.9627,38:4.956,40:4.949}
def im_at(N,D):
    w=sla.eigs(liou(N,D,Om,V),k=4,sigma=-0.894,return_eigenvectors=False,tol=1e-13); w=w[np.argsort(abs(w+0.894))][:2]
    return abs(w[0].imag), w
def job(N):
    Ds=np.linspace(c0[N]-0.006,c0[N]+0.006,121); ims=np.array([im_at(N,D)[0] for D in Ds]); m=ims>1e-8
    if not m.any(): return [N,None,None,float(ims.max())]
    return [N,float(Ds[m].min()),float(Ds[m].max()),float(ims.max())]
if __name__=="__main__":
    with Pool(7) as P: res=P.map(job,sorted(c0))
    for r in res: print(r,flush=True)
    json.dump(res,open('../data/rwin2.json','w'))
