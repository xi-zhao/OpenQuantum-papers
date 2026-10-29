import numpy as np, scipy.sparse.linalg as sla, json, sys
from scipy.optimize import minimize_scalar
from ex import liou
from multiprocessing import Pool
V=10.
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
def pair(N,D,Om,t):
    w=sla.eigs(liou(N,D,Om,V),k=4,sigma=t,return_eigenvectors=False,tol=1e-13)
    return np.sort(w[np.argsort(abs(w-t))][:2].real)
def job(N):
    D0=EP[N][1]; t=EP[N][3]
    res=minimize_scalar(lambda D: np.diff(pair(N,D,2.0,t))[0],bracket=(D0-0.005,D0,D0+0.005),tol=1e-11)
    return [N,res.x,res.fun/2]
if __name__=="__main__":
    with Pool(8) as P: out=P.map(job,[24,28,32,36,40,44,48,52,56])
    json.dump(out,open('../data/gN.json','w'))
    for r in out: print(r,flush=True)
