import numpy as np, scipy.sparse.linalg as sla
from ex import liou
V=10.;Om=2.0
for N in [32,36]:
    for D in np.linspace(4.90,5.05,16):
        w=sla.eigs(liou(N,D,Om,V),k=6,sigma=-0.89,return_eigenvectors=False,tol=1e-12); w=w[np.argsort(abs(w+0.89))][:3]
        print(N,round(D,3),np.round(w,5),flush=True)
